struct Bilin_interp {
	Int m,n;
	const MatDoub &y;
	Linear_interp x1terp, x2terp;
		
	Bilin_interp(VecDoub_I &x1v. VecDoub_I &x2v, MatDoub_I &ym)
		: m(x1v.size()), n(x2v.size()), y(ym),
		x1terp(x1v,x1v), x2terp(x2v,x2v) {}
	/*
	Construct dummy 1-dim interpolations ofr their locate and hunt functions.
	*/
	Doub interp(Doub x1p, Doub x2p) {
		Int i,j;
		Doub yy, t, u;
		i = x1terp.cor ? x1terp.hunt(x1p) : x1terp.locate(x1p);
		j = x2terp.cor ? x2terp.hunt(x2p) : x2terp.locate(x2p);
		//find the grid square.
		t = (x1p-x1terp.xx[i]) / (x1terp.xx[i+1]-x1terp.xx[i]);
		u = (x2p-x2terp.xx[j]) / (x2terp.xx[j+1]-x2terp.xx[j]);
		yy = (1.-t)*(1.-u)*y[i][j] + t*(1.-u)*y[i+1][j]
			+(1.-t)*u*y[i][j+1] + t*u*y[i+1][j+1];
		return yy;
	}
};
//\TODO data type MatDoub and VecDoub need to be defined. From Numerical Recipes
