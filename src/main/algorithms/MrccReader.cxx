/* Copyright 2021 cc4s.org
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <algorithms/MrccReader.hpp>

#include <tcc/Tcc.hpp>
#include <Complex.hpp>
#include <MathFunctions.hpp>
#include <Log.hpp>
#include <Exception.hpp>
#include <Cc4s.hpp>
#include <extern/Lapack.hpp>
#include <TensorSet.hpp>


using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(MrccReader)

Ptr<MapNode> MrccReader::run(const Ptr<MapNode> &arguments) {
  Ptr<MapNode> result;
  // multiplex calls to template methods
  if (Cc4s::dryRun) {
    using TE = DefaultDryTensorEngine;
    result = run<TE>(arguments);
  } else {
    using TE = DefaultTensorEngine;
    result = run<TE>(arguments);
  }
  ASSERT_LOCATION(
    result, "unsupported tensor type as 'operator'",
    arguments->sourceLocation
  );
  return result;
}

template <typename TE>
Ptr<MapNode> MrccReader::run(
  const Ptr<MapNode> &arguments
) {
  auto result(New<MapNode>(SOURCE_LOCATION));

/*
cpr(i) is just the largest element of MO(i,mu) for all mu

*/

  // map of fock matrix elements in AO basis
  auto _fockAO  = fileRead<Data2>("FOCK");
  // map of MO-Coefficients
  auto _fullMos = fileRead<Data2>("MOCOEF");
  // index map needed for sparse AO CoulombVertex
  auto intpos   = fileRead<Data3>("TEDAT");
  // We have to change the indices from

//  if (Cc4s::dryRun) return result;
  auto rank(Cc4s::world->getRank());

  int ii(0);
  for (auto &i: intpos) {
    i.index1--;
    i.index2--;
  }
  for (auto &d: _fockAO) {
    d.index1--;
    d.index2--;
  }
  for (auto &d: _fullMos) {
    d.index1--;
    d.index2--;
  }


  // we need to have access to these numbers
  size_t nocc = 4;
  size_t nvir = 64;
  size_t ncore = 1;
  size_t nbasis = nocc + nvir;
  size_t nbf = nbasis + ncore; //28
  size_t ndfbasis = 152;
  std::vector<Real<>> fockAO(nbf*nbf), fockMO(nbasis*nbasis);
  std::vector<Real<>> mos(nbasis*nbf), tmp(nbasis*nbf);
  std::vector<Real<>> eps(nbasis), epsp(nvir), epsh(nocc);
  std::vector<Real<>> mos_occ(nbf*nocc), mos_vir(nbf*nvir);

  for (auto &f: _fockAO) {
    fockAO[f.index1 + nbf * f.index2] = f.value;
    fockAO[f.index2 + nbf * f.index1] = f.value;
  }

  size_t jjj(0);
  for (auto &m: _fullMos) {
    if (m.index2 < ncore) continue;
    mos[m.index1 + nbf * (m.index2 - ncore)] = m.value;
    if (m.index2 < ncore+nocc) {
      mos_occ[m.index1 + nbf * (m.index2 - ncore)] = m.value;
    } else {
      mos_vir[m.index1 + nbf * (m.index2 - ncore - nocc)] = m.value;
    }
  }



  cc4s::Real<64> alpha = 1.0, beta = 0.0;
  dsymm_("l",
         "l",
         (int const *) &nbf,
         (int const *) &nbasis,
         &alpha,
         fockAO.data(),
         (int const *) &nbf,
         mos.data(),
         (int const *) &nbf,
         &beta,
         tmp.data(),
         (int const *) &nbf);

  dgemm_("t",
         "n",
         (int const *) &nbasis,
         (int const *) &nbasis,
         (int const *) &nbf,
         &alpha,
         mos.data(),
         (int const *) &nbf,
         tmp.data(),
         (int const *) &nbf,
         &beta,
         fockMO.data(),
         (int const *) &nbasis);

  size_t p(0), h(0);
  for (size_t i(0); i < nbasis; i++) {
    eps[i] = fockMO[i + nbasis*i];
    if (i <  nocc) epsh[h++] = fockMO[i + nbasis*i];
    if (i >= nocc) epsp[p++] = fockMO[i + nbasis*i];
  }

  auto halfVertexh = halfTransform(mos_occ, intpos, nocc, nbf, ndfbasis);

  auto Mos_vir(Tcc<TE>::template tensor<Real<>>({nbf, nvir}, "Mos_vir") );
  auto Gamma_igm(Tcc<TE>::template tensor<Real<>>({nocc, ndfbasis, nbf}, "halfVertexh"));
  auto Epsp(Tcc<TE>::template tensor<Real<>>({nvir}, "epsp") );
  auto Epsh(Tcc<TE>::template tensor<Real<>>({nocc}, "epsh") );

  size_t el = (rank) ? 0 : halfVertexh.size();
  // Note: we use the 'large' index map for all tensors. watch out!
  std::vector<size_t> idx(el);
  std::iota(idx.begin(), idx.end(), 0);
  Gamma_igm->write(el, idx.data(), halfVertexh.data());

  el = (rank) ? 0 : mos_vir.size();
  Mos_vir->write(el, idx.data(), mos_vir.data());
  el = (rank) ? 0 : nvir;
  Epsp->write(el, idx.data(), epsp.data());
  el = (rank) ? 0 : nocc;
  Epsh->write(el, idx.data(), epsh.data());

  auto _vrtx = Gamma_igm->readAll();

  auto gamma_ab = eval_mp2(Gamma_igm, Epsp, Epsh, Mos_vir, true);

  int lwork(3*nvir), info;
  std::vector<Real<>> w(nvir), work(lwork);

  dsyev_("V",
         "L",
         (int const *) &nvir,
         gamma_ab.data(),
         (int const *) &nvir,
         w.data(),
         work.data(),
         (int const *) &lwork,
         (int const *) &info);
  if (info != 0) throw "problem diagonalization (1), naturalOrbitals\n";

  size_t fnoFactor(arguments->template getValue<Natural<>>("fnoFactor"));
  size_t nFno = std::min(fnoFactor*nocc, nvir);

  for (size_t a(0); a < nvir - nFno; a++)
  for (size_t b(0); b < nvir; b++) {
    gamma_ab[b + a * nvir] = 0.0;
  }


  std::vector<Real<>> fock_ab(nvir*nvir), scratch(nvir*nvir), rot_ab(nFno*nvir);
  //construct new fock matrix using the truncated density matrix
  for (size_t c(0); c < nvir; c++)
  for (size_t b(0); b < nvir; b++)
  for (size_t a(0); a < nvir; a++) {
    fock_ab[a + b * nvir] += gamma_ab[c + a * nvir] * gamma_ab[c + b * nvir] * epsp[c];
  }

  dsyev_("V",
         "L",
         (int const *) &nvir,
         fock_ab.data(),
         (int const *) &nvir,
         w.data(),
         work.data(),
         (int const *) &lwork,
         (int const *) &info);
  if (info != 0) throw "problem diagonalization (4), naturalOrbitals\n";


  el = (rank) ? 0 : nFno;
  auto Epsf( Tcc<TE>::template tensor<Real<>>({nFno}, "Epsf"));
  auto Mos_vir_fno( Tcc<TE>::template tensor<Real<>>({nbf, nFno}, "Mos_vir_fno"));
  Epsf->write(el, idx.data(), &w[nvir-nFno]);

  for (size_t b(0); b < nvir; b++)
  for (size_t a(0); a < nvir; a++)
  for (size_t c(0); c < nvir; c++) {
    scratch[a + b * nvir] += gamma_ab[a + c * nvir] * fock_ab[c + b * nvir];
  }
  for (size_t a(0); a < nvir; a++) {
    size_t c(0);
    for (size_t b(nvir-nFno); b < nvir; b++) {
      rot_ab[a+ nvir*c++ ] = scratch[a+b*nvir];
    }
  }

  // ~C[mu,b] = R[a,b] * C[mu,a]
  std::vector<Real<>> mos_vir_fno(nFno*nbf);
  for (size_t b(0); b < nFno; b++)
  for (size_t m(0); m < nbf; m++)
  for (size_t a(0); a < nvir; a++) {
    mos_vir_fno[m + b * nbf] += mos_vir[m + a * nbf] * rot_ab[a + b * nvir];
  }

  el = (rank) ? 0 : nFno*nbf;
  Mos_vir_fno->write(el, idx.data(), mos_vir_fno.data());

  eval_mp2(Gamma_igm, Epsf, Epsh, Mos_vir_fno, false);

  auto halfVertexp = halfTransform(mos_vir_fno, intpos, nFno, nbf, ndfbasis);

  auto Gamma_agm(Tcc<TE>::template tensor<Real<>>({nFno, ndfbasis, nbf}, "halfVertexp"));
  el = (rank) ? 0 : halfVertexp.size();
  idx.resize(el);
  std::iota(idx.begin(), idx.end(), 0);
  Gamma_agm->write(el, idx.data(), halfVertexp.data());


  auto Mos_occ(Tcc<TE>::template tensor<Real<>>({nbf, nocc}, "Mos_occ"));
  el = (rank) ? 0 : nocc*nbf;
  Mos_occ->write(el, idx.data(), mos_occ.data());


  auto Gamma_gij(Tcc<TE>::template tensor<Complex<>>({ndfbasis, nocc, nocc}, "gamma_gij"));
  auto Gamma_gia(Tcc<TE>::template tensor<Complex<>>({ndfbasis, nocc, nFno}, "gamma_gia"));
  auto Gamma_gai(Tcc<TE>::template tensor<Complex<>>({ndfbasis, nFno, nocc}, "gamma_gai"));
  auto Gamma_gab(Tcc<TE>::template tensor<Complex<>>({ndfbasis, nFno, nFno}, "gamma_gab"));

  COMPILE(
    (*Gamma_gij)["Gij"] <<= map<Complex<>>(toComplex<Real<>>, (*Gamma_igm)["iGm"] * (*Mos_occ)["mj"]),
    (*Gamma_gia)["Gia"] <<= map<Complex<>>(toComplex<Real<>>, (*Gamma_igm)["iGm"] * (*Mos_vir_fno)["ma"]),
    (*Gamma_gai)["Gai"] <<= (*Gamma_gia)["Gia"],
    (*Gamma_gab)["Gab"] <<= map<Complex<>>(toComplex<Real<>>, (*Gamma_agm)["aGm"] * (*Mos_vir_fno)["mb"])
  )->execute();

  std::shared_ptr<TensorSet<Complex<>,TE>> vertex = std::make_shared<TensorSet<Complex<>,TE>>();
  vertex->get("hh") = Gamma_gij;
  vertex->get("hp") = Gamma_gia;
  vertex->get("ph") = Gamma_gai;
  vertex->get("pp") = Gamma_gab;

  std::shared_ptr<TensorSet<Real<>,TE>> energies = std::make_shared<TensorSet<Real<>,TE>>();
  energies->get("h") = Epsh;
  energies->get("p") = Epsf;


  result->setPtr("coulombVertex", vertex);
  result->setPtr("eigenEnergies", energies);



  return result;

}


template <typename T>
std::vector<T> MrccReader::fileRead(const std::string& filename, size_t nint) {
  std::ifstream iffile(filename, std::ios::binary);
  if (!iffile) {
    // TODO we have to cc4s assert this
    //throw std::runtime_error("Failed to open file: " + filename);
  }

  // Skip Fortran 4-byte record marker at the beginning
  iffile.seekg(4, std::ios::beg);

  // if nint is not given we read the first entry of the file
  if (!nint) {
    iffile.read(reinterpret_cast<char*>(&nint), sizeof(size_t));
  }
  // Allocate vector and read data
  std::vector<T> data(nint);
  iffile.read(reinterpret_cast<char*>(data.data()), sizeof(T) * nint);

  // Ensure we read the expected amount
  if (iffile.gcount() != static_cast<std::streamsize>(sizeof(T) * nint)) {
    // TODO we have to cc4s assert this
    //throw std::runtime_error("Error reading binary file: " + filename);
  }

  return data;
}

std::vector<Real<>> MrccReader::halfTransform(
  const std::vector<Real<>> &mos,
  const std::vector<Data3>  &intpos,
  const size_t norb,
  const size_t nbf,
  const size_t ndfbasis
) {

  std::ifstream iffile("TEINT", std::ios::binary);
  iffile.seekg(0, std::ios::beg);
  std::vector<Real<>> halfVertex(norb*nbf*ndfbasis);
  std::vector<Real<>> motr(norb*nbf);
  std::vector<Real<>> kvec(ndfbasis);

  for (size_t i(0); i < norb; i++)
  for (size_t m(0); m < nbf; m++) {
    motr[i + m*norb] = mos[m + i* nbf];
  }
  for (auto &ipos: intpos) {
    std::fill(kvec.begin(), kvec.end(), 0);
    auto m(ipos.index1);
    auto n(ipos.index2);
    auto nn(ipos.index3);
    //read n elements from the file
    std::vector<Data1> line(nn);


    iffile.read(reinterpret_cast<char*>(line.data()), sizeof(Data1) * nn);
    for (auto &k: line) {
      // fotran to C index shift
      kvec[k.index-1] = k.value;
    }
    for (size_t g(0); g < ndfbasis; g++)
    for (size_t i(0); i < norb; i++) {
      // \Gamma["iG\nu"] = C["i\mu"] * kvec["G"]
      // \Gamma["iG\mu"] = C["i\nu"] * kvec["G"]
      halfVertex[i + g* norb + n*norb*ndfbasis] += motr[i + m*norb] * kvec[g];
      if (m == n) continue;
      halfVertex[i + g* norb + m*norb*ndfbasis] += motr[i + n*norb] * kvec[g];
    }
  }

//  for (size_t n(0); n < dims.nbf; n++)
//  for (size_t g(0); g < dims.ndfbasis; g++) {
//    auto idx = g* dims.nocc + n*dims.nocc*dims.ndfbasis;
//    if (halfVertex[idx] * halfVertex[idx] > 1e-15)
//      OUT() << n << " " << g << " | " << halfVertex[idx] << std::endl;
//  }

  return halfVertex;
}

template <typename TE>
std::vector<Real<>> MrccReader::eval_mp2(Ptr<Tensor<Real<>,TE>> Gamma_igm,
                                         Ptr<Tensor<Real<>,TE>> Epsp,
                                         Ptr<Tensor<Real<>,TE>> Epsh,
                                         Ptr<Tensor<Real<>,TE>> Mos_vir,
                                         bool eval_gamma_ab) {

  auto Gamma_ab(Tcc<TE>::template tensor<Real<>>("gamma_ab"));
  auto Energy(Tcc<TE>::template tensor<Real<>>("energy"));

  auto _lens(Gamma_igm->getLens());
  auto ndfbasis(_lens[1]);
  auto nocc(_lens[0]);
  auto nvir(Epsp->getLen(0));

  auto Vpphh(Tcc<TE>::template tensor<Real<>>("Vpphh"));
  auto Dpphh(Tcc<TE>::template tensor<Real<>>({nvir, nvir, nocc, nocc}, "Dpphh"));
  auto Tpphh(Tcc<TE>::template tensor<Real<>>("Tpphh"));
  auto Gamma_gia(Tcc<TE>::template tensor<Real<>>({ndfbasis, nocc, nvir}, "Vertex_ph"));
  COMPILE(
    (*Gamma_gia)["Gia"] <<= (*Gamma_igm)["iGm"] * (*Mos_vir)["ma"],
    (*Vpphh)["abij"] <<= (*Gamma_gia)["Gia"] * (*Gamma_gia)["Gjb"],
    (*Dpphh)["abij"]  -= (*Epsp)["a"],
    (*Dpphh)["abij"]  -= (*Epsp)["b"],
    (*Dpphh)["abij"]  += (*Epsh)["i"],
    (*Dpphh)["abij"]  += (*Epsh)["j"]
  )->execute();
  Dpphh->sum(1.0, "abij", "abij", [] (Real<> delta) { return 1.0 / delta; });
  COMPILE(
    //(*energy)[""] <<= 2.0 * (*Vpphh)["abij"] * (*Vpphh)["abij"]
    (*Energy)[""] <<= 2.0 * (*Vpphh)["abij"] * (*Vpphh)["abij"] * (*Dpphh)["abij"],
    (*Energy)[""]  -= 1.0 * (*Vpphh)["abij"] * (*Vpphh)["baij"] * (*Dpphh)["abij"],
    // Now we evaluate the Mp2 density matrix //
    (*Tpphh)["abij"]  <<= 2.0 * (*Vpphh)["abij"] * (*Dpphh)["abij"],
    (*Tpphh)["abij"]   -= 0.5 * (*Tpphh)["baij"],
    (eval_gamma_ab) ? (
      (*Gamma_ab)["ab"] <<= (*Tpphh)["caij"] * (*Vpphh)["cbij"] * (*Dpphh)["cbij"],
      Tcc<TE>::sequence()
    ) : (
      Tcc<TE>::sequence()
    )
  )->execute();

  auto E(Energy->read());
  std::cout.precision(10);
  OUT() <<   E << std::endl;

  if (!eval_gamma_ab) return std::vector<Real<>>();

  auto el = (!Cc4s::world->getRank()) ? nvir*nvir : 0;
  std::vector<size_t> idx(el);
  std::iota(idx.begin(), idx.end(), 0);
  std::vector<Real<>> gamma_ab(el);
  Gamma_ab->read(el, idx.data(), gamma_ab.data());

  return gamma_ab;
}


