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
#include <lapack/LapackWrapper.hpp>
#include <TensorSet.hpp>
#include <string>


using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(MrccReader)

Ptr<MapNode> MrccReader::run(const Ptr<MapNode> &arguments) {
  Ptr<MapNode> result;
  // multiplex calls to template methods
  if (Cc4s::dryRun) {
    using TE = DefaultDryTensorEngine;
    result = run<TE, DryMode::DRY>(arguments);
  } else {
    using TE = DefaultTensorEngine;
    result = run<TE, DryMode::WET>(arguments);
  }
  ASSERT_LOCATION(
    result, "unsupported tensor type as 'operator'",
    arguments->sourceLocation
  );
  return result;
}




// Get Γab specialization for Dry Tensors (DefaultDryTensorEngine, DryMode::DRY)
template <>
std::vector<Real<>>
MrccReader::get_gamma_ab<DefaultDryTensorEngine, DryMode::DRY>(bool fnoApprox,
                                                               ___tensor_wrapper<DefaultDryTensorEngine> Gamma_igm,
                                                               ___tensor_wrapper<DefaultDryTensorEngine> Epsp,
                                                               ___tensor_wrapper<DefaultDryTensorEngine> Epsh,
                                                               ___tensor_wrapper<DefaultDryTensorEngine> Mos_vir) {
  if (fnoApprox) {
    using DryTensor = Tensor<Real<>, DefaultDryTensorEngine>;
    return eval_mp2_approx_dry<DefaultDryTensorEngine>(Gamma_igm,
                                                       Epsp,
                                                       Epsh,
                                                       Mos_vir);
  } else {
    return eval_mp2(Gamma_igm, Epsp, Epsh, Mos_vir, true);
  }
}

// Get Γab specialization for Wet Tensors (DefaultTensorEngine, DryMode::WET)
template <>
std::vector<Real<>> MrccReader::get_gamma_ab<DefaultTensorEngine, DryMode::WET>(bool fnoApprox,
                                                                                ___tensor_wrapper<DefaultTensorEngine> Gamma_igm,
                                                                                ___tensor_wrapper<DefaultTensorEngine> Epsp,
                                                                                ___tensor_wrapper<DefaultTensorEngine> Epsh,
                                                                                ___tensor_wrapper<DefaultTensorEngine> Mos_vir) {
  if (fnoApprox) {
    using WetTensor = Tensor<Real<>, DefaultTensorEngine>;
    return eval_mp2_approx<DefaultTensorEngine>(Gamma_igm,
                                                Epsp,
                                                Epsh,
                                                Mos_vir);
  } else {
    return eval_mp2(Gamma_igm, Epsp, Epsh, Mos_vir, true);
  }
}

template <typename TE, DryMode dry_mode>
Ptr<MapNode> MrccReader::run(
  const Ptr<MapNode> &arguments
) {
  auto result(New<MapNode>(SOURCE_LOCATION));


  size_t nocc, nvir, ncore, ndfbasis, nbasis, nbf;
  if (arguments->isGiven("mrccStdOut")) {
    auto mrccOut(arguments->template getValue<std::string>("mrccStdOut"));
    auto systemDimensions = extractSystemDimensions(mrccOut);
    ncore = systemDimensions.core_electrons /2;
    nocc = systemDimensions.electrons / 2 - ncore;
    ndfbasis = systemDimensions.correlation_fitting;
    nbf = systemDimensions.conventionalBasis;
    nbasis = nbf - ncore;
    nvir = nbasis - nocc;
  } else {
    nocc = arguments->template getValue<Natural<>>("nocc");
    nvir = arguments->template getValue<Natural<>>("nvir");
    ncore = arguments->template getValue<Natural<>>("ncore");
    ndfbasis = arguments->template getValue<Natural<>>("ndfbasis");
    nbasis = nocc + nvir;
    nbf = nbasis + ncore;
  }
  size_t fnoFactor(arguments->template getValue<Natural<>>("fnoFactor", 0));
  bool fnoApprox(arguments->template getValue<int>("fnoApprox", 0) == 1);
  if (fnoApprox) {
    OUT() << "Using approximated Frozen Natural Orbitals (FNOs)\n";
  }

  size_t nFno = std::min(fnoFactor*nocc, nvir);
  if (fnoFactor == 0) nFno = nvir;

  std::vector<Real<>> mos_occ(nbf*nocc), mos_vir_fno(nFno*nbf);
  std::vector<Real<>> mos_vir, epsp, epsh;

  // index map needed for sparse AO CoulombVertex
  auto intpos   = fileRead<Data3>("TEDAT");
  for (auto &i: intpos) {
    i.index1--;
    i.index2--;
  }

  size_t rank(Cc4s::world->getRank());
  auto np(Cc4s::getProcessesCount());
  OUT() << "Mrcc full vertex: nocc: " << nocc << " nvir: " << nvir
        << " nbf: " << nbf << " ndf: " << ndfbasis << std::endl;
  LOG() << "Mrcc full vertex: nocc: " << nocc << " nvir: " << nvir
        << " nbf: " << nbf << " ndf: " << ndfbasis << std::endl;

  if (!rank) {

    // map of fock matrix elements in AO basis
    auto _fockAO  = fileRead<Data2>("FOCK");
    // map of MO-Coefficients
    auto _fullMos = fileRead<Data2>("MOCOEF");

    // We have to change the indices from fortran to C
    for (auto &d: _fockAO) {
      d.index1--;
      d.index2--;
    }
    for (auto &d: _fullMos) {
      d.index1--;
      d.index2--;
    }


    // we need to have access to these numbers
    std::vector<Real<>> fockAO(nbf*nbf), fockMO(nbasis*nbasis);
    std::vector<Real<>> mos(nbasis*nbf), tmp(nbasis*nbf);
    std::vector<Real<>> eps(nbasis);
    epsp.resize(nvir);
    epsh.resize(nocc);
    mos_vir.resize(nvir*nbf);
    for (auto &f: _fockAO) {
      fockAO[f.index1 + nbf * f.index2] = f.value;
      fockAO[f.index2 + nbf * f.index1] = f.value;
    }

    for (auto &m: _fullMos) {
      if ((size_t) m.index2 < ncore) continue;
      mos[m.index1 + nbf * (m.index2 - ncore)] = m.value;
      if ((size_t) m.index2 < ncore+nocc) {
        mos_occ[m.index1 + nbf * (m.index2 - ncore)] = m.value;
      } else {
        mos_vir[m.index1 + nbf * (m.index2 - ncore - nocc)] = m.value;
      }
    }

    matmul(fockAO, mos, tmp, nbf, nbasis, nbf);
    matmul(mos, tmp, fockMO, nbasis, nbasis, nbf, true, false);

    size_t p(0), h(0);
    for (size_t i(0); i < nbasis; i++) {
      eps[i] = fockMO[i + nbasis*i];
      if (i <  nocc) epsh[h++] = fockMO[i + nbasis*i];
      if (i >= nocc) epsp[p++] = fockMO[i + nbasis*i];
    }
  }

  // we have to distribute mos_occ over the different ranks
  // right now only rank 0 is used
  //auto halfVertexh = halfTransform(mos_occ, intpos, nocc, nbf, ndfbasis)
  LOG() << "Reading and transforming three center integrals" << std::endl;
  DryMemory::allocate(ndfbasis*nbf*nocc, sizeof(double), "halfVertexh", SOURCE_LOCATION);
  std::vector<Real<>> halfVertexh;
  if (Cc4s::dryRun) {
    std::ifstream iffile("TEINT", std::ios::binary);
    ASSERT_LOCATION(
      iffile, "File TEINT not present",
      SOURCE_LOCATION
    );
    // check the file length of TEINT
    int64_t te_size(0);
    for (auto ipos: intpos) {
      auto nn(ipos.index3);
      te_size += nn * sizeof(Data1);
    }
    LOG() << "FILE SIZE OF TEINT IS " << te_size <<  " BYTES" << std::endl;
  } else {
    MPI_Bcast(mos_occ.data(), nbf*nocc, MPI_DOUBLE, 0, Cc4s::world->getComm());
    halfVertexh = halfTransform(mos_occ, intpos, nocc, nbf, ndfbasis);
  }
  auto Mos_vir(Tcc<TE>::template tensor<Real<>>({nbf, nvir}, "Mos_vir") );
  auto Gamma_igm(Tcc<TE>::template tensor<Real<>>({nocc, ndfbasis, nbf}, "halfVertexh"));
  auto Epsp(Tcc<TE>::template tensor<Real<>>({nvir}, "epsp") );
  auto Epsh(Tcc<TE>::template tensor<Real<>>({nocc}, "epsh") );

  size_t el;
  DryMemory::allocate(ndfbasis*nbf*nocc, sizeof(size_t), "halfVertexh_idx", SOURCE_LOCATION);
  std::vector<size_t> idx;
  if (!Cc4s::dryRun) {
    el = halfVertexh.size();
    idx.resize(el);
    auto nbf_loc = nbf / np + (rank < nbf % np);
    for (size_t m(0); m < nbf_loc; m++)
    for (size_t g(0); g < ndfbasis; g++)
    for (size_t i(0); i < nocc; i++) {
      auto m_global = m * np + rank;
      idx[i + g*nocc + m * nocc * ndfbasis] = i + g*nocc + m_global * nocc * ndfbasis;
    }
    Gamma_igm->write(el, idx.data(), halfVertexh.data());
  }
  // clear std::vector
  halfVertexh = std::vector<Real<>>();
  idx = std::vector<size_t>();
  DryMemory::free(ndfbasis*nbf*nocc, sizeof(double), "halfVertexh");
  DryMemory::free(ndfbasis*nbf*nocc, sizeof(double), "halfVertexh_idx");

  el = (rank) ? 0 : mos_vir.size();
  idx.resize(el);
  std::iota(idx.begin(), idx.end(), 0);
  Mos_vir->write(el, idx.data(), mos_vir.data());
  el = (rank) ? 0 : nvir;
  Epsp->write(el, idx.data(), epsp.data());
  el = (rank) ? 0 : nocc;
  Epsh->write(el, idx.data(), epsh.data());


  DryMemory::allocate(np*nvir*nvir, sizeof(double), "gamma_ab_glb", SOURCE_LOCATION);
  auto gamma_ab = get_gamma_ab< TE, dry_mode >(fnoApprox, Gamma_igm, Epsp, Epsh, Mos_vir);


  std::vector<Real<>> w(nvir);
  if (!Cc4s::dryRun) {
    if (!rank) {
      LOG() << "eval FNO: Diagonalize density matrix" << std::endl;
      compute_eigen(gamma_ab, nvir, w);
      // zero out the elements beyond nFno
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

      LOG() << "eval FNO: Diagonalize Fock matrix" << std::endl;
      compute_eigen(fock_ab, nvir, w);

      for (size_t b(0); b < nvir; b++)
      for (size_t a(0); a < nvir; a++)
      for (size_t c(0); c < nvir; c++) {
        scratch[a + b * nvir] += gamma_ab[a + c * nvir] * fock_ab[c + b * nvir];
      }
      //matmul(scratch, gamma_ab, fock_ab, nvir, nvir, nvir);

      for (size_t a(0); a < nvir; a++) {
        size_t c(0);
        for (size_t b(nvir-nFno); b < nvir; b++) {
          rot_ab[a+ nvir*c++ ] = scratch[a+b*nvir];
        }
      }


      // ~C[mu,b] = R[a,b] * C[mu,a]
      for (size_t b(0); b < nFno; b++)
      for (size_t m(0); m < nbf; m++)
      for (size_t a(0); a < nvir; a++) {
        mos_vir_fno[m + b * nbf] += mos_vir[m + a * nbf] * rot_ab[a + b * nvir];
      }
      //matmul(mos_vir_fno, mos_vir, rot_ab, nbf, nFno, nvir, true, true);

    }
    // we have to distribute to all ranks for the AO->MO transform
    LOG() << "Distribute AO->MO transform among all ranks" << std::endl;
    MPI_Bcast(mos_vir_fno.data(), nbf*nFno, MPI_DOUBLE, 0, Cc4s::world->getComm());
  }

  gamma_ab = std::vector<Real<>>();
  DryMemory::free(np*nvir*nvir, sizeof(double), "gamma_ab_glb");

  el = (rank) ? 0 : nFno;

  auto Epsf( Tcc<TE>::template tensor<Real<>>({nFno}, "Epsf"));
  auto Mos_vir_fno( Tcc<TE>::template tensor<Real<>>({nbf, nFno}, "Mos_vir_fno"));
  Epsf->write(el, idx.data(), &w[nvir-nFno]);

  el = (rank) ? 0 : nFno*nbf;
  idx.resize(el);
  std::iota(idx.begin(), idx.end(), 0);
  Mos_vir_fno->write(el, idx.data(), mos_vir_fno.data());


  if (!fnoApprox) {
    eval_mp2(Gamma_igm, Epsf, Epsh, Mos_vir_fno, false);
  }
  std::vector<Real<>> halfVertexp;

  DryMemory::allocate(nFno*ndfbasis*nbf, sizeof(double), "halfVertexp", SOURCE_LOCATION);

  LOG() << "Half transform virtual-virtual Vertex" << std::endl;
  if (!Cc4s::dryRun) halfVertexp = halfTransform(mos_vir_fno, intpos, nFno, nbf, ndfbasis);

  auto Gamma_agm(Tcc<TE>::template tensor<Real<>>({nFno, ndfbasis, nbf}, "halfVertexp"));

  DryMemory::allocate(nFno*ndfbasis*nbf, sizeof(size_t), "halfVertexp_idx", SOURCE_LOCATION);
  if (!Cc4s::dryRun) {
    el = halfVertexp.size();
    idx.resize(el);
    auto nbf_loc = nbf / np + (rank< nbf % np);
    for (size_t m(0); m < nbf_loc; m++)
    for (size_t g(0); g < ndfbasis; g++)
    for (size_t a(0); a < nFno; a++) {
      auto m_global = m * np + rank;
      idx[a + g*nFno + m * nFno * ndfbasis] = a + g*nFno + m_global * nFno * ndfbasis;
    }
  }
  Gamma_agm->write(el, idx.data(), halfVertexp.data());

  DryMemory::free(nFno*ndfbasis*nbf, sizeof(double), "halfVertexp");
  halfVertexp = std::vector<Real<>>();
  DryMemory::free(nFno*ndfbasis*nbf, sizeof(size_t), "halfVertexp_idx");
  idx = std::vector<size_t>();


  auto Mos_occ(Tcc<TE>::template tensor<Real<>>({nbf, nocc}, "Mos_occ"));
  el = (rank) ? 0 : nocc*nbf;
  idx.resize(el);
  std::iota(idx.begin(), idx.end(), 0);

  Mos_occ->write(el, idx.data(), mos_occ.data());


  auto Gamma_gij(Tcc<TE>::template tensor<Real<>>({ndfbasis, nocc, nocc}, "gamma_gij"));
  auto Gamma_gia(Tcc<TE>::template tensor<Real<>>({ndfbasis, nocc, nFno}, "gamma_gia"));
  auto Gamma_gai(Tcc<TE>::template tensor<Real<>>({ndfbasis, nFno, nocc}, "gamma_gai"));
  auto Gamma_gab(Tcc<TE>::template tensor<Real<>>({ndfbasis, nFno, nFno}, "gamma_gab"));

  COMPILE(
    (*Gamma_gij)["Gij"] <<= (*Gamma_igm)["iGm"] * (*Mos_occ)["mj"],
    (*Gamma_gia)["Gia"] <<= (*Gamma_igm)["iGm"] * (*Mos_vir_fno)["ma"],
    (*Gamma_gai)["Gai"] <<= (*Gamma_gia)["Gia"],
    (*Gamma_gab)["Gab"] <<= (*Gamma_agm)["aGm"] * (*Mos_vir_fno)["mb"]
  )->execute();

  std::shared_ptr<TensorSet<Real<>,TE>> vertex = std::make_shared<TensorSet<Real<>,TE>>();
  vertex->get("hh") = Gamma_gij;
  vertex->get("hp") = Gamma_gia;
  vertex->get("ph") = Gamma_gai;
  vertex->get("pp") = Gamma_gab;

  std::shared_ptr<TensorSet<Real<>,TE>> energies = std::make_shared<TensorSet<Real<>,TE>>();
  energies->get("h") = Epsh;
  energies->get("p") = Epsf;

  result->setPtr("coulombVertex", vertex);
  result->setPtr("eigenEnergies", energies);


  context["occ"] = nocc;
  context["vir"] = nFno;
  context["NF"]  = ndfbasis;
  context["type"] = "Real64";



  return result;

}


template <typename T>
std::vector<T> MrccReader::fileRead(const std::string& filename, size_t nint) {
  std::ifstream iffile(filename, std::ios::binary);
  ASSERT_LOCATION(
    iffile, "File " + filename + " not present",
    SOURCE_LOCATION
  );

  // Skip Fortran 4-byte record marker at the beginning
  iffile.seekg(4, std::ios::beg);

  // if nint is not given we read the first entry of the file
  if (!nint) {
    iffile.read(reinterpret_cast<char*>(&nint), sizeof(size_t));
  }
  // Allocate vector and read data
  std::vector<T> data(nint);
  iffile.read(reinterpret_cast<char*>(data.data()), sizeof(T) * nint);

  ASSERT_LOCATION(
    iffile.gcount() == static_cast<std::streamsize>(sizeof(T) * nint),
    "File " + filename + " not present",
    SOURCE_LOCATION
  );

  return data;
}

std::pair<bool,bool> isForMe(int64_t m, int64_t n) {
  auto rank(Cc4s::world->getRank());
  auto np(Cc4s::world->getProcesses());
  return { m % np == rank, n % np == rank};
}

std::vector<Real<>> MrccReader::halfTransform(
  const std::vector<Real<>> &mos,
  const std::vector<Data3>  &intpos,
  const size_t norb,
  const size_t nbf,
  const size_t ndfbasis
) {
  MPI_File fh;
  MPI_Status status;
  MPI_File_open(MPI_COMM_WORLD, "TEINT", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
//  std::ifstream iffile("TEINT", std::ios::binary);
//  iffile.seekg(0, std::ios::beg);

  size_t rank(Cc4s::world->getRank());
  auto np(Cc4s::world->getProcesses());
  auto nbf_loc = nbf / np + (rank< nbf % np);
  std::vector<Real<>> halfVertex(norb*nbf_loc*ndfbasis);
  //std::vector<Real<>> halfVertex(norb*nbf*ndfbasis);
  std::vector<Real<>> motr(norb*nbf);
  std::vector<Real<>> kvec(ndfbasis);

  for (size_t i(0); i < norb; i++)
  for (size_t m(0); m < nbf; m++) {
    motr[i + m*norb] = mos[m + i* nbf];
  }
  MPI_Offset offset = 0;
  for (auto &ipos: intpos) {
    std::fill(kvec.begin(), kvec.end(), 0);
    auto m(ipos.index1);
    auto n(ipos.index2);
    auto nn(ipos.index3);
    auto b(isForMe(m,n));
    // is this rank does not any data we still have to update the offset
    if (!b.first && !b.second) {
      offset += nn*sizeof(Data1);
      continue;
    }
    //read n elements from the file
    std::vector<Data1> line(nn);
    MPI_File_read_at(fh, offset, line.data(), nn * sizeof(Data1), MPI_BYTE, &status);
    offset += nn*sizeof(Data1);
    //iffile.read(reinterpret_cast<char*>(line.data()), sizeof(Data1) * nn);
    for (auto &k: line) {
      // fotran to C index shift
      kvec[k.index-1] = k.value;
    }
    for (size_t g(0); g < ndfbasis; g++)
    for (size_t i(0); i < norb; i++) {
      // \Gamma["iG\nu"] = C["i\mu"] * kvec["G"]
      // \Gamma["iG\mu"] = C["i\nu"] * kvec["G"]
/* this is the non-MPI version
      halfVertex[i + g* norb + n*norb*ndfbasis] += motr[i + m*norb] * kvec[g];
      if (m == n) continue;
      halfVertex[i + g* norb + m*norb*ndfbasis] += motr[i + n*norb] * kvec[g];
*/
      size_t m_loc = m / np, n_loc = n /np;
      if (b.first) {
        halfVertex[i + g* norb + m_loc*norb*ndfbasis] += motr[i + n*norb] * kvec[g];
      }
      if (b.second && m != n) {
        halfVertex[i + g* norb + n_loc*norb*ndfbasis] += motr[i + m*norb] * kvec[g];
      }
    }
  }

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
  //auto Tpphh(Tcc<TE>::template tensor<Real<>>("Tpphh"));
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
    // D = V/\Delta
    (*Dpphh)["abij"] <<= (*Vpphh)["abij"] * (*Dpphh)["abij"],
    (*Energy)[""] <<= 2.0 * (*Vpphh)["abij"] * (*Dpphh)["abij"],
    (*Energy)[""]  -= 1.0 * (*Vpphh)["abij"] * (*Dpphh)["baij"],
    (eval_gamma_ab) ? (
      (*Gamma_ab)["ab"] <<= 2.0 * (*Dpphh)["caij"] * (*Dpphh)["cbij"],
      (*Gamma_ab)["ab"]  -= 1.0 * (*Dpphh)["caij"] * (*Dpphh)["cbji"],
      Tcc<TE>::sequence()
    ) : (
      Tcc<TE>::sequence()
    )
/*
    (*Energy)[""] <<= 2.0 * (*Vpphh)["abij"] * (*Vpphh)["abij"] * (*Dpphh)["abij"],
    (*Energy)[""]  -= 1.0 * (*Vpphh)["abij"] * (*Vpphh)["baij"] * (*Dpphh)["abij"],
    // Now we evaluate the Mp2 density matrix //
    (eval_gamma_ab) ? (
      (*Tpphh)["abij"]  <<= 2.0 * (*Vpphh)["abij"] * (*Dpphh)["abij"],
      (*Tpphh)["abij"]   -= 0.5 * (*Tpphh)["baij"],
      (*Gamma_ab)["ab"] <<= (*Tpphh)["caij"] * (*Vpphh)["cbij"] * (*Dpphh)["cbij"],
      Tcc<TE>::sequence()
    ) : (
      Tcc<TE>::sequence()
    )
*/
  )->execute();

  if (!eval_gamma_ab) return std::vector<Real<>>();

  auto E(Energy->read());
  std::cout.precision(10);
  OUT() << "Mp2 correlation energy(full basis): " <<   E << std::endl;
  auto el = (!Cc4s::world->getRank()) ? nvir*nvir : 0;
  std::vector<size_t> idx(el);
  std::iota(idx.begin(), idx.end(), 0);
  std::vector<Real<>> gamma_ab(el);
  Gamma_ab->read(el, idx.data(), gamma_ab.data());

  return gamma_ab;
}

template <typename TE>
std::vector<Real<>> MrccReader::eval_mp2_approx_ref(Ptr<Tensor<Real<>,TE>> Gamma_igm,
                                                    Ptr<Tensor<Real<>,TE>> Epsp,
                                                    Ptr<Tensor<Real<>,TE>> Epsh,
                                                    Ptr<Tensor<Real<>,TE>> Mos_vir) {



  auto Gamma_ab(Tcc<TE>::template tensor<Real<>>("gamma_ab"));

  auto _lens(Gamma_igm->getLens());
  auto ndfbasis(_lens[1]);
  auto nocc(_lens[0]);
  auto nvir(Epsp->getLen(0));

  auto Vpphh(Tcc<TE>::template tensor<Real<>>("Vpphh"));
  auto Dpphh(Tcc<TE>::template tensor<Real<>>({nvir, nvir, nocc, nocc}, "Dpphh"));
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
    (*Dpphh)["abij"]  <<= (*Vpphh)["abij"] * (*Dpphh)["abij"],
    (*Gamma_ab)["ab"] <<= (*Dpphh)["cbii"] * (*Dpphh)["caii"]
  )->execute();

  auto el = (!Cc4s::world->getRank()) ? nvir*nvir : 0;
  std::vector<size_t> idx(el);
  std::iota(idx.begin(), idx.end(), 0);
  std::vector<Real<>> gamma_ab(el);
  Gamma_ab->read(el, idx.data(), gamma_ab.data());

  return gamma_ab;
}

template <typename TE>
std::vector<Real<>> MrccReader::eval_mp2_approx(Ptr<Tensor<Real<>,TE>> Gamma_igm,
                                                Ptr<Tensor<Real<>,TE>> Epsp,
                                                Ptr<Tensor<Real<>,TE>> Epsh,
                                                Ptr<Tensor<Real<>,TE>> Mos_vir) {

  auto Gamma_ab(Tcc<TE>::template tensor<Real<>>("gamma_ab"));

  auto _lens(Gamma_igm->getLens());
  auto ndfbasis(_lens[1]);
  auto nocc(_lens[0]);
  auto nvir(Epsp->getLen(0));

  auto epsp = Epsp->readAll();
  auto epsh = Epsh->readAll();
  auto Vpphh(Tcc<TE>::template tensor<Real<>>("Vpphh"));
  auto Dpphh(Tcc<TE>::template tensor<Real<>>({nvir, nvir, nocc, nocc}, "Dpphh"));
  auto Gamma_gia(Tcc<TE>::template tensor<Real<>>({ndfbasis, nocc, nvir}, "Vertex_ph"));
  COMPILE(
    (*Gamma_gia)["Gia"] <<= (*Gamma_igm)["iGm"] * (*Mos_vir)["ma"]
  )->execute();
  std::vector<double> gamma_ab_glb;
  // this is a CTF::Tensor now
  auto _Gamma_gia = Gamma_gia->evaluate()->getMachineTensor()->tensor;
  // we want to evaluate:
  // \gamma_ab = \sum_i [ \sum_c ( \sum_g \Gamma[gci] * \Gamma[gbi] / \Delta_bc^ii)
  //                           * ( \sum_g \Gamma[gci] * \Gamma[gai] / \Delta_ac^ii) ]
  // for each i: buf_ga = \Gamma[gci]
  // buf_pp = \sum_g \Gamma[gci] * \Gamma[gbi]
  //
  // 1.) we distribute Gamma over the last index and calculate different i on
  //     different processors (some processors can hold more than one occupied index
  // 2.) we evaluate a local result of \gamma_ab and accumulate the result in the end


  int np(Cc4s::world->getProcesses()), rank(Cc4s::world->getRank());
  std::vector<int64_t> lens{static_cast<int64_t>(ndfbasis),
                            static_cast<int64_t>(nvir),
                            static_cast<int64_t>(nocc)};
  std::vector<int> sym(3, NS);
  std::vector<int> plens(3, 1);
  plens[2] = np;
  CTF::Partition part(3, plens.data());
  CTF::Tensor<> G_gai(3, lens.data(), sym.data(), *_Gamma_gia.wrld, "Gai", part["Gai"]);

  G_gai["Gai"] = _Gamma_gia["Gia"];

  std::vector<double> gamma_ab(nvir*nvir, 0.0);

  size_t localOcc = nocc / np + size_t( rank < nocc % np);
  size_t b_size(nvir*ndfbasis);
  std::vector<double> Dabii(nvir*nvir), buf_ga(b_size), buf_pp(nvir*nvir);
  for (size_t ii(0UL); ii < localOcc; ii++) {
    size_t occ(ii*np + rank);
    double epsi(epsh[occ]);
    for (size_t v(0UL); v < nvir; v++) {
      for (size_t w(0UL); w < nvir; w++) {
        Dabii[v * nvir + w] = 2*epsi - epsp[w] - epsp[v];
      }
    }
    double a(1.0), b(0.0);
    int inc(1), v(nvir), g(ndfbasis);
    size_t offset = ii * b_size * sizeof(double);
    std::memcpy(buf_ga.data(), G_gai.data + offset, b_size*sizeof(double));
    dgemm_("T", "N", &v, &v, &g, &a, buf_ga.data(), &g, buf_ga.data(), &g, &b, buf_pp.data(), &v);
    for (size_t v(0UL); v < buf_pp.size(); v++) buf_pp[v] /= Dabii[v];
    dgemm_("N", "T", &v, &v, &v, &a, buf_pp.data(), &v, buf_pp.data(), &v, &a, gamma_ab.data(), &v);
  }
  MPI_Barrier(Cc4s::world->getComm());
  // accumulate the result to rank 0 who carries the full density matrix
  gamma_ab_glb.resize(nvir*nvir, 0.0);
  MPI_Reduce(
    gamma_ab.data(),
    gamma_ab_glb.data(),
    nvir * nvir,
    MPI_DOUBLE,
    MPI_SUM,
    0,
    Cc4s::world->getComm()
  );


  return gamma_ab_glb;
}

template <typename TE>
std::vector<Real<>> MrccReader::eval_mp2_approx_dry(Ptr<Tensor<Real<>,TE>> Gamma_igm,
                                                    Ptr<Tensor<Real<>,TE>> Epsp,
                                                    Ptr<Tensor<Real<>,TE>> Epsh,
                                                    Ptr<Tensor<Real<>,TE>> Mos_vir) {

  auto _lens(Gamma_igm->getLens());
  auto ndfbasis(_lens[1]);
  auto nocc(_lens[0]);
  auto nvir(Epsp->getLen(0));

  auto epsp = Epsp->readAll();
  auto epsh = Epsh->readAll();
  auto Vpphh(Tcc<TE>::template tensor<Real<>>("Vpphh"));
  auto Dpphh(Tcc<TE>::template tensor<Real<>>({nvir, nvir, nocc, nocc}, "Dpphh"));
  auto Gamma_gia(Tcc<TE>::template tensor<Real<>>({ndfbasis, nocc, nvir}, "Vertex_ph"));
  COMPILE(
    (*Gamma_gia)["Gia"] <<= (*Gamma_igm)["iGm"] * (*Mos_vir)["ma"]
  )->execute();

  size_t np(Cc4s::getProcessesCount());

  DryMemory::allocate(ndfbasis*nocc*nvir, sizeof(double), "Gamma_gia_partitioned", SOURCE_LOCATION);
  DryMemory::allocate(np*nvir*nvir, sizeof(double), "gamma_ab", SOURCE_LOCATION);
  DryMemory::allocate(np*nvir*nvir, sizeof(double), "Dabii", SOURCE_LOCATION);
  DryMemory::allocate(np*nvir*nvir, sizeof(double), "buf_pp", SOURCE_LOCATION);
  DryMemory::allocate(np*nvir*ndfbasis, sizeof(double), "buf_ga", SOURCE_LOCATION);

  DryMemory::free(ndfbasis*nocc*nvir, sizeof(double), "Gamma_gia_partitioned");
  DryMemory::free(np*nvir*nvir, sizeof(double), "gamma_ab");
  DryMemory::free(np*nvir*nvir, sizeof(double), "Dabii");
  DryMemory::free(np*nvir*nvir, sizeof(double), "buf_pp");
  DryMemory::free(np*nvir*ndfbasis, sizeof(double), "buf_ga");


  return std::vector<double>();
}


// Extract a number from a string like "Total number of basis functions:  xxx"
int extractNumberFromLine(const std::string& line) {
    size_t colonPos = line.find(':');
    if (colonPos == std::string::npos) return 0;

    std::string numberPart = line.substr(colonPos + 1);
    // Remove all whitespace
    numberPart.erase(
        std::remove_if(numberPart.begin(), numberPart.end(),
                      [](unsigned char c) { return std::isspace(c); }),
        numberPart.end());

    // Convert to integer
    return std::stoi(numberPart);
}

MrccReader::SystemDimensions MrccReader::extractSystemDimensions(const std::string filename) {
    std::ifstream infile(filename);
    ASSERT_LOCATION(
      infile, "Mrcc Stdout not present. File " + filename ,
      SOURCE_LOCATION
    );

    LOG() << "Start reading the Mrcc StdOut." << std::endl;
    std::string line;
    MrccReader::SystemDimensions result = {0, 0, 0, 0, -1};

    // Track which section we're in
    enum Section { NONE, CONVENTIONAL, SCF_FITTING, CORRELATION_FITTING };
    Section currentSection = NONE;

    // Track which block we're in
    int blockCount = 0;
    //bool inTargetBlock = useFirstBlock; // If true, use first block; otherwise wait for second
    bool useFirstBlock = false;
    bool inTargetBlock = useFirstBlock; // If true, use first block; otherwise wait for second

    while (std::getline(infile, line)) {
        // Detect section changes
        if (line.find("Conventional basis set:") != std::string::npos) {
            currentSection = CONVENTIONAL;
            if (currentSection == CONVENTIONAL && blockCount == 0) {
                blockCount++;
                inTargetBlock = useFirstBlock;
            } else if (currentSection == CONVENTIONAL && blockCount == 1) {
                blockCount++;
                inTargetBlock = !useFirstBlock;
            }
        } else if (line.find("Fitting basis set for the SCF calculation:") != std::string::npos) {
            currentSection = SCF_FITTING;
        } else if (line.find("Fitting basis set for the correlation calculation:") != std::string::npos) {
            currentSection = CORRELATION_FITTING;
        }

        // Check if this is the line with basis function count
        if (inTargetBlock && line.find("Total number of basis functions:") != std::string::npos) {
            int value = extractNumberFromLine(line);

            switch (currentSection) {
                case CONVENTIONAL:
                    result.conventionalBasis = value;
                    break;
                case SCF_FITTING:
                    result.scf_fitting = value;
                    break;
                case CORRELATION_FITTING:
                    result.correlation_fitting = value;
                    break;
                default:
                    break;
            }
        }
        if (line.find(" Number of electrons:") != std::string::npos) {
          result.electrons = extractNumberFromLine(line);
        }
        if (line.find("umber of core electrons:") != std::string::npos) {
          result.core_electrons = extractNumberFromLine(line);
        }


        // If we've filled all values and we're in the target block, we can stop
        if (inTargetBlock &&
            result.conventionalBasis != 0 &&
            result.correlation_fitting != 0 &&
            result.electrons != 0 &&
            result.core_electrons >= 0) {
            break;
        }

    }
    ASSERT_LOCATION(
      result.conventionalBasis &&
      result.correlation_fitting &&
      result.electrons,
      "Mrcc Stdout not read sucessfully",
      SOURCE_LOCATION
    );
    return result;
}

