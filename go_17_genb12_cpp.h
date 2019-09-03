void GEN_BANDES_12::SGUA2::Debug(const char * lib) {
	cout << lib << " gua2 \t" << i_81 << " cols " << col1 + 1 << col2 + 1
		<< " digs " << dig1 + 1 << dig2 + 1 
		<< " digs 0"<<oct<<digs<<dec<< endl;
}

void GEN_BANDES_12::InitialSockets2Setup() {//load permanent data
	for (int i = 0; i < 81; i++) {// initial socket 2
		SGUA2 & w= tsgua2[i];
		w.i_81 = i;
		register int i9 = i/9;
		int tpcol[3][2] = { {1,2},{0,2},{0,1} };
		int rdcol = i9 % 3,dcol=3*(i9/3),*p=tpcol[rdcol];
		w.i9 = i9;
		w.col1 = dcol + p[0];
		w.col2 = dcol + p[1];
		int rd1 = i - 9 * w.i9;//relative d1;d2		
		w.id1 = rd1 / 3 + 3 * w.col1;
		w.id2 = rd1 % 3 + 3 * w.col2;

	}
}
void GEN_BANDES_12::InitialSockets3Setup() {//load permanent data
	for (int i = 0; i < 81; i++) {// initial socket 3
		SGUA3 & w = tsgua3[i];
		int minir(i / 27);
		w.i_81 = i;
		w.imini = minir;
		w.col1 = 3 * minir;// minirow first column in gangster
		int dp = i - 27 * minir,// the 9 perms of the gangster minirow
			dg = 9 * minir;// gangster 27 start
		 // pointers to gangster digits
		w.id1 = dp / 9;
		dp -= 9 * w.id1;
		w.id1 += dg;
		dg += 3;// next gangster column inthe minirow
		w.id2 = dp / 3;
		dp -= 3 * w.id2;
		w.id2 += dg;
		dg += 3;// last gangster column in the minirow
		w.id3 = dp + dg;
	}
}

void GEN_BANDES_12::SecondSockets2Setup() {
	ntua2 = 0; nactive2 = 0;
	Build_CheckUAs_Subsets_Table();
	for (i81 = 0; i81 < 81; i81++) {// initial socket 2
		SGUA2 & w = tsgua2[i81];
		w.dig1 = gang27[w.id1];
		w.dig2 = gang27[w.id2];
		w.digs = (1 << w.dig1) | (1 << w.dig2);
		if (!gchk.IsGua(i81)) continue;
		// build revised gangster
		memcpy(w.gangcols, gangb12, sizeof gangb12);
		w.gangcols[w.col1] ^= w.digs;
		w.gangcols[w.col2] ^= w.digs;
		// find guas of the socket
		zh2b_g.InitGangster(gangb12, w.gangcols);
		zh2b5_g.sizef5 = GUALIMSIZE;
		zh2b5_g.modevalid = 0;
		w.nua = 0;
		ptua2 = w.tua;
		nua2 = 0;

		//================== GUA collector 2 bands 
		GuaCollect(w.digs );
		for (int i = 0; i < 84; i++) {// find UAs 3 digits
			int fl = floors_3d[i];
			if ((fl & w.digs) == w.digs)	GuaCollect(fl);
		}
		for (int i = 0; i < 126; i++) {// find UAs 4 digits
			int fl = floors_4d[i];
			if ((fl & w.digs) == w.digs) GuaCollect(fl);
		}
		for (int i = 0; i < 126; i++) {// find UAs 5digits
			int fl = 0x1ff ^ floors_4d[i];
			if ((fl & w.digs) == w.digs) GuaCollect(fl);
		}
		SecondSockets2MoreUAs();
		if (nua2) {
			tactive2[nactive2++]=i81;
			ntua2 += nua2;
		}
		w.nua = nua2;
	}
}
void GEN_BANDES_12::SecondSockets3Setup() {
	ntua3 = 0; nactive3 = 0;
	for (int i81 = 0; i81 < 81; i81++) {// initial socket 2
		SGUA3 &w = tsgua3[i81];
		w.dig1 = gang27[w.id1];
		w.dig2 = gang27[w.id2];
		w.dig3 = gang27[w.id3];
		if (!gchk.IsGua3(i81))continue;
		// Setup the perms for gangsters in minirow
		int bita = 1 << w.dig1, bitb = 1 << w.dig2, bitc = 1 << w.dig3,
			digs = bita | bitb | bitc;
		int triplet_perms[2][3];

		int * p = triplet_perms[0];// for per abc -> bca
		p[0] = bita | bitb; p[1] = bitb | bitc; p[2] = bitc | bita;

		p = triplet_perms[1];// for per abc -> cab
		p[0] = bita | bitc; p[1] = bitb | bita; p[2] = bitc | bitb;
		ptua2 = w.tua;
		nua2 = 0;

		// build revised gangster
		//w.gangcols[w.col2] ^= w.digs;
		int tp3f[2][3] = { {1,2,0},{2,0,1} };// perms no valid digit
		for (int ip = 0; ip < 2; ip++) {
			// build revised gangster
			int rgangcols[9];// revised gangster
			memcpy(rgangcols, gangb12, sizeof gangcols);
			p = triplet_perms[ip];
			int c1 = w.col1, c2 = c1 + 1, c3 = c1 + 2;
			rgangcols[c1] ^= p[0];
			rgangcols[c2] ^= p[1];
			rgangcols[c3] ^= p[2];

			// find guas of the socket
			zh2b_g.InitGangster(gangb12, rgangcols);
			zh2b5_g.sizef5 = UALIMSIZE - 1;//; - 3;
			zh2b5_g.modevalid = 0;
			//w.nua_start = ntua3;

			//================== GUA collector 2 bands 
			GuaCollect(digs);// find UAs 3 digits
			for (int i = 0; i < 126; i++) {// find UAs 4 digits
				int fl = floors_4d[i];
				if ((fl& digs) == digs)	GuaCollect(fl);
			}
			for (int i = 0; i < 126; i++) {// find UAs 5digits
				int fl = 0x1ff ^ floors_4d[i];
				if ((fl & digs) == digs) GuaCollect(fl);
			}
		}
		if (nua2) {
			tactive3[nactive3++] = i81;
			ntua3 += nua2;
		}
		w.nua = nua2;// store the final count
	}
}

void ZH2B::Init_2digits_banda(BF64  cellsbf) {
	int txcells[50], nxcells = 0;
	BitsInTable64(txcells,nxcells, cellsbf.bf.u64);
	memset(zh2b_g.Digit_cell_Assigned_init, 0, sizeof zh2b_g.Digit_cell_Assigned_init);
	for (int i = 0; i < nxcells; i++) {
		int xcell=txcells[i],cell=From_128_To_81[xcell],
		digit = zh2b_g.puz0[cell];
		Assign(digit, cell, xcell);
		zh2b_g.Digit_cell_Assigned_init[digit].Set(xcell);
	}
	for (int i = 0; i < 9; i++)  FD[i] &= cells_unsolved |
		zh2b_g.Digit_cell_Assigned_init[i];
}

/*
extract from bands 1 2 and UAs bands1+2 all uas (uint64_t)
having a chance to be a subset of a gua
Note : this is to save time. A UA having a subset is not an error.
Note 2 : All digits of a subset are in the UA to check  
*/
void GEN_BANDES_12::Build_CheckUAs_Subsets_Table() {
	nua_for_check = 0;
	for (uint32_t i = 0; i < genuasb12.nuab1b2; i++) {// uAs bands
		register uint64_t R = genuasb12.tuab1b2[i]& BIT_SET_2X;
		uint64_t cc = _popcnt64(R);
		if (cc > 10) continue;// nearly no chance to be subset
		tua_for_check[nua_for_check++] = R| (cc<<59);
	}
	for (uint32_t i = 0; i < genuasb12.nua; i++){
		register uint64_t R = genuasb12.tua[i];
		if ((R >> 59) > 10) continue;// nearly no chance to be subset
		tua_for_check[nua_for_check++] = R;
	}
	// and build the digit pattern of the uas in tua_for_check
	for (uint32_t i = 0; i < nua_for_check; i++) {
		uint32_t digs = 0,res;
		register uint64_t R = tua_for_check[i]&BIT_SET_2X;
		while (bitscanforward64(res, R)) {
			R ^= (uint64_t)1 << res;
			int cell = From_128_To_81[res],digit=gchk.grid0[cell];
			digs |= 1 << digit;
		}
		uadigs_for_check[i] = digs;
		if (0) {
			cout << Char2Xout(tua_for_check[i]) << " ua for check digs 0"
				<< oct << digs << dec <<" ndigs "<<_popcnt32(digs)<< endl;
		}
	}

}
void GEN_BANDES_12::Build_tuacheck(int fl) {
	nua_check=0;
	int ndigs = 0x1ff ^ fl;
	for (uint32_t i = 0; i < nua_for_check; i++) {
		uint32_t digs = uadigs_for_check[i];
		if (digs & ndigs) continue;
		tuacheck[nua_check++]= tua_for_check[i];
	}
}
int GEN_BANDES_12::Have_tuacheck_subset(  ) {// ua 2x27  + 5 bit length
	uint64_t * t = tuacheck;
	register uint64_t ua2x = ua & BIT_SET_2X;
	for (uint32_t iua = 0; iua < nua_check; iua++) {
		register uint64_t R = t[iua];
		R &= BIT_SET_2X;
		if ((R&ua2x) == R)		return 1;// we have a subset
	}
	return 0;
}

void GEN_BANDES_12::SecondSockets2MoreUAs() {
	SGUA2 s = tsgua2[i81];
	int rx0[2], bx0[2], dx3[2];
	rx0[0] = zh2b_g.fd_sols[0][s.dig1].SolRow(s.col2);
	dx3[0] = gchk.grid0[9 * rx0[0] + s.col1];
	rx0[1] = zh2b_g.fd_sols[0][s.dig2].SolRow(s.col1);
	dx3[1] = gchk.grid0[9 * rx0[1] + s.col2];
	bx0[0] = rx0[0] / 3;
	bx0[1] = rx0[1] / 3;

	if (bx0[0] == bx0[1]) return; // want digits in 2 bands
	if (dx3[0] == dx3[1]) return; // ua 6 cells
	for (int idig = 0; idig < 2; idig++) {
		int dx[2], cx[2],row,dig3;
		if (idig) {
			dx[0] = s.dig1; dx[1] = s.dig2;
			cx[0] = s.col1; cx[1] = s.col2;
			row = rx0[0]; dig3 = dx3[0];
		}
		else {
			dx[0] = s.dig2; dx[1] = s.dig1;
			cx[0] = s.col2; cx[1] = s.col1;
			row = rx0[1]; dig3 = dx3[1];
		}
		band=row/3;
		// dig3 must be possible in 
		if (!((1 << dig3) & s.gangcols[cx[1]])) {
			continue;
		}
		// now fill band bx[0] except the 2 cells
		zh2b[0].Init_gang();
		BF64 cells_to_fill;
		cells_to_fill.bf.u64=BIT_SET_2X;// all cells
		cells_to_fill.bf.u32[1-band] = 0;// nothing in the other band
		c1 = row * 9 + s.col1;
		c2 = row * 9 + s.col2;
		if (row >= 3) { c1 += 5; c2 += 5; }// must be 2X mode
		cells_to_fill.bf.u64 ^= (uint64_t)1 << c1;
		cells_to_fill.bf.u64 ^= (uint64_t)1 << c2;
		zh2b[0].Init_2digits_banda(cells_to_fill);
		zh2b[0].FullUpdate();
		if(diagmore)zh2b[0].ImageCandidats();
		// at this point,we have only one band with a gangster 
		//not equal to the start and a partial GUA 2cells in the other band
		//we  prepare ZHONE to find partial uas in one band
		uint32_t tua[300]; // supply a tua to ZHONE
		zh1b_g.tua = tua;
		for (int i = 0; i < 9; i++) {// load the gangster pm
			zh1b_g.fd_sols[1][i] = zh2b[0].FD[i].bf.u32[1 - band];
		}
		STD_B416 * myb = &gchk.bands_abc[0];
		if(!band) myb= &gchk.bands_abc[1]; // to catch the solution
		memcpy(zh1b_g.fd_sols[0], myb->fd_sols[0], sizeof zh1b_g.fd_sols[0]);
		zh1b_g.band0 = myb->band0;
		digp = s.digs | (1 << dig3);// must be in the floor

		//==============  now find more uas 6/7 digits
		ngua6_7 = 6;
		wua0 = ((uint64_t)1 << c1) | ((uint64_t)1 << c2);
		for (int i6 = 0; i6 < 84; i6++) {
			int fl3 = floors_3d[i6];// , fl6 = 0x1ff ^ fl3;
			if (fl3&digp) continue;// digits must be in the multi floors
			floors = 0x1ff ^ fl3;
			GuaCollectMore();
		}
		ngua6_7 = 7;
		for (int i7 = 0; i7 < 36; i7++) {
			int fl2 = floors_2d[i7];// , fl7 = 0x1ff ^ fl2;
			if (fl2&digp) continue;// digits must be in the multi floors
			floors = 0x1ff ^ fl2;
			GuaCollectMore();
		}
	}
}
void GEN_BANDES_12::GuaCollectMore() {
	//if (diagmore)if (floors == 0663) { diagmore = 2; zh1b_g.diag = 1; }
	//else { diagmore = 1; zh1b_g.diag = 0; }
	zh1b_g.nua = 0;
	zhone[0].Start_nFloors(floors);
	zhone[0].InitGuess();
	if (diagmore > 1) zhone[0].ImageCandidats();
	else zh1b_g.diag = 0;
	if (ngua6_7 == 6)zhone[0].Guess6();
	else zhone[0].Guess7();
	if (zh1b_g.nua) {
		Build_tuacheck(floors);
		SGUA2 s = tsgua2[i81];
		for (uint32_t i = 0; i < zh1b_g.nua; i++) {
			uint32_t w = zh1b_g.tua[i]&BIT_SET_27,	cc=_popcnt32(w),cell;
			//cout << Char27out(w) << endl;;
			if (cc < 8) continue;
			if (cc < 12 && ngua6_7 == 7) continue;
			if (cc > GUALIMSIZE - 2) continue;
			int digpw = digp;// assigned outside the band
			register uint32_t R = w;
			while (bitscanforward(cell, R)) {
				R ^= 1 << cell;
				int dig = zh1b_g.band0[cell]; // digit changed
				digpw |= 1 << dig;
			}
			if (_popcnt32(digpw) != ngua6_7) continue;
			//cout << "valid ua to add to the table" << endl;
			ua = wua0;
			ua |= (uint64_t)w << 32 * (1 - band);
			if (Have_tuacheck_subset()) {
				//cout << Char2Xout(ua) << " ua has subset" << endl;
				continue;
			}
			ua |= (uint64_t)(cc + 2) << 59;
			if (nua2 >= SIZETGUA)nua2 = SIZETGUA - 1; // guess it will be a smaller
			int ir = 		AddUA64(ptua2, nua2, ua);
			if(ir && diagmore)		cout << Char2Xout(ua) << " wua added cc=" << cc + 2 <<" nua2="<<nua2<< endl;

		}
	}
	// try to add the ua to the file
}
void GEN_BANDES_12::GuaCollect(int fl,int diag) {//use revised gangster
	uint64_t solved_cells = zh2b5_g.FindUAsInit(fl, 1);
	if (!solved_cells) return;// one digit solved true
	if (diag) cout << "return from zh2b5_g.FindUAsInit(fl, 1);" << endl;
	zh2b5_g.CollectUas5();// collect uas for this set of floors
	if (diag) cout << "zh2b5_g.CollectUas5();" << endl;
	if (!zh2b5_g.nuaf5) return;
	Build_tuacheck( fl);

	// check subsets and add to main table
	for (uint32_t i = 0; i < zh2b5_g.nuaf5; i++) {
		ua = zh2b5_g.tuaf5[i].bf.u64;
		//genuasb12.ua = ua;
		if (diag)cout << Char2Xout(ua) << " to add if new2 nua=" << nua2 << endl;
		if (Have_tuacheck_subset()) continue;// superset of a previous ua
		//		if (genuasb12.CheckOld()) continue;// superset of a previous ua
		uint64_t cc = _popcnt64(ua&BIT_SET_2X);
		ua |= cc << 59;
		//if (genuasb12.CheckMain(ua)) continue;// superset of a previous ua
		// see why, seems to be of no effect
		// protect against table limit
		if (nua2 >= SIZETGUA)nua2 = SIZETGUA - 1; // guess it will be a smaller
		int ir = AddUA64(ptua2, nua2,ua);
		if (ir) {
			if (diag)cout << Char2Xout(genuasb12.ua) <<" ir="<<ir << " added nua2="<<nua2 << endl;
		}
	}
}

void GEN_BANDES_12::ValidInitGang() {
	for (int i = 0; i < 9; i++) {// init columns status
		cold[i] = 0x1ff;
		for (int j = 0; j < 6; j++)	cold[i] ^= 1 << gchk.grid0[i + 9 * j];
		gangb12[i] = 0x1ff ^ cold[i];
	}
	memcpy(gangcols, cold, sizeof gangcols);
}
void GEN_BANDES_12::BuildGang9x3() {
	gang27 = gang[0];
	for (int i = 0; i < 9; i++) {// 9 cols to build out of gangcols
		int istack = C_box[i];
		int * d = gang[i], c = gangcols[i];
		uint32_t bit;
		bitscanforward(bit, c);
		c ^= (1 << bit);
		d[0] = bit;
		gang_digits_cols[bit][istack] = i;
		bitscanforward(bit, c);
		c ^= (1 << bit);
		d[1] = bit;
		gang_digits_cols[bit][istack] = i;
		bitscanforward(bit, c);
		d[2] = bit;
		gang_digits_cols[bit][istack] = i;
	}

}
void GEN_BANDES_12::DebugGang() {

	cout << "gang 12 gang 3< gang digits" << endl;
	for (int i = 0; i < 9; i++) {
		cout << oct << "0" << gangb12[i] << "\t";
	}
	cout <<endl;
	for (int i = 0; i < 9; i++) {
		cout << oct << "0" << gangcols[i] << "\t";
	}
	cout << endl;
	for (int i = 0; i < 9; i++) {
		cout << oct << "0" << (gangcols[i]| gangb12[i]) << "\t";
	}
	cout << endl;
	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 3; j++) cout <<  gang[i][j]  ;
		cout <<   "\t";
	}
	cout << endl;
}

int GEN_BANDES_12::DebugFreshUA(uint64_t ua) {
	// purely debugging code, a fresh UA <= limsize must have >5 digits
	uint32_t cell, digits = 0;
	while (bitscanforward64(cell, ua)) {
		ua ^= (uint64_t)1 << cell;
		digits |= 1 << gchk.grid0[From_128_To_81[cell]];
	}
	if (_popcnt32(digits) >= 6) return 0;
	cout << "bug more ua bands 1+2 not enough digits " << endl;
	return 1;
}
