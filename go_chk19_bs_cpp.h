void GCHK::Start(STD_B416 * bandx, int * tp, int * tsort, int mode) {
	{ // catch the three band in the right order and mode
		for (int ib = 0; ib < 3; ib++) {
			int iperm = tp[ib], iband = tsort[iperm] & 7;
			bands_abc[ib] = bandx[iband];
			memcpy(&zsol[27 * ib], bands_abc[ib].band, 27);
		}
		mode_equal = mode; // 0-3 to tell where = clues have to be considered
		// set up ordered solution
		zsol[81] = 0;
		for (int i = 0; i < 81; i++)grid0[i] = zsol[i] - '1';
		cout << zsol << " traité" << endl;
		genb12.ValidInitGang();// set up gang bandc and gang bandsab
		genb12.BuildGang9x3();// setup ordered gang digis band C
		{// Init band3 guas process 
			memset(&guas, 0, sizeof guas);
			int * band0 = &grid0[54];
			for (int i = 0; i < 9; i++) {
				minirows_bf[i] = 0;
				int * p = &band0[3 * i];
				for (int j = 0; j < 3; j++)
					minirows_bf[i] |= 1 << p[j];
			}
		}
	}
	DoExpandBand();// expand band A all valid bands 3-6 clues
	//DebugExpand();
	//=========================== collect UAs  
	zh1b_g.modegua = 0;//must be to activate filter in UAs b12 more
	if (genuasb12.Initgen()) return;
	zh1b_g.modegua = 1;//must be to kill  filter in GUAs 6_7 more
	genb12.SecondSockets2Setup();// collect GUA2s 
	genb12.SecondSockets3Setup();// collect GUA3s 
	if(0)GodebugInit(0);
	GoBanda();
}
void GCHK::DoExpandBand() {
	uint32_t * tua = bands_abc[0].tua, nua = bands_abc[0].nua;
	struct SPB3 {// spots to find band 3 minimum valid solutions
		// ====================== constant after initialization
		int  possible_cells, all_previous_cells, active_cells, iuab3;
		//uint64_t cursol;
	}spb3[8], *s3, *sn3;
	s3 = spb3;
	s3->all_previous_cells = 0;
	s3->active_cells = BIT_SET_27;// all cells active
	s3->iuab3 = 0; // copy the start table
	s3->possible_cells = tua[0];
	n3_6 = 0;
	nxindex3 = 0;
	int tcells[10];
	//____________________  here start the search
next:
	uint64_t ispot = s3 - spb3;
	uint32_t cell;
	if (!bitscanforward(cell, s3->possible_cells))goto back;
	{// apply cell in bitfields
		register int bit = (uint64_t)1 << cell;
		tcells[ispot] = cell;
		s3->possible_cells ^= bit;// clear bit
		register int filter = s3->all_previous_cells | bit,
			ac = s3->active_cells ^ bit;
		sn3 = s3 + 1;
		sn3->all_previous_cells = filter;
		if (ispot == 2) 	xindex3[nxindex3++].Open(filter, n3_6, tcells);
		sn3->active_cells = s3->active_cells = ac;
		// nextspot:take the next available ua to loop
		for (int i = s3->iuab3 + 1; i < (int)nua; i++) {
			if (tua[i] & filter)continue;
			if (ispot >= 5) 	goto next;//passing the limit
			sn3->iuab3 = i;
			sn3->possible_cells = tua[i] & ac;
			s3 = sn3; // switch to next spot
			goto next;
		}

	}
	// no more ua
	if (ispot == 5)	x_expand_3_6[n3_6++].Add3(tcells[3], tcells[4], tcells[5]);
	else {// if below 5 loop  for redundant clues
		int t[32], nt = 0;
		register int ac = s3->active_cells;
		while (bitscanforward(cell, ac)) {// put active cells in table
			register int bit = 1 << cell;
			ac ^= bit;
			t[nt++] = cell;
		}
		if (ispot == 4) {// valid 5 clues
			x_expand_3_6[n3_6++].Add2(tcells[3], tcells[4]);
			for (int i6 = 0; i6 < nt; i6++)
				x_expand_3_6[n3_6++].Add3(tcells[3], tcells[4], t[i6]);
		}
		//cout << "nt="<<nt<<" ispot="<<ispot << endl;
		else if (ispot == 3) { // valid 4 clues
			x_expand_3_6[n3_6++].Add1(tcells[3]);
			for (int i5 = 0; i5 < nt; i5++) {
				x_expand_3_6[n3_6++].Add2(tcells[3], t[i5]);
				for (int i6 = i5+1; i6 < nt; i6++)
					x_expand_3_6[n3_6++].Add3(tcells[3], t[i5], t[i6]);
			}
		}
		else if (ispot == 2) { // valid 3 clues
			x_expand_3_6[n3_6++].Add0();
			for (int i4 = 0; i4 < nt; i4++) {
				x_expand_3_6[n3_6++].Add1(t[i4]);
				tcells[3] = t[i4];
				for (int i5 = i4 + 1; i5 < nt; i5++) {
					x_expand_3_6[n3_6++].Add2(t[i4], t[i5]);
					for (int i6 = i5 + 1; i6 < nt; i6++)
						x_expand_3_6[n3_6++].Add3(t[i4],  t[i5], t[i6]);
				}
			}
		}
		else if (ispot == 1) { // valid 2 clues
			for (int i3 = 0; i3 < nt; i3++) {
				tcells[2] = t[i3];
				int filter = (1 << tcells[0]) | (1 << tcells[1]) | (1 << t[i3]);
				xindex3[nxindex3++].Open(filter, n3_6, tcells);
				x_expand_3_6[n3_6++].Add0();
				for (int i4 = i3 + 1; i4 < nt; i4++) {
					x_expand_3_6[n3_6++].Add1(t[i4]);
					for (int i5 = i4 + 1; i5 < nt; i5++) {
						x_expand_3_6[n3_6++].Add2(t[i4], t[i5]);
						for (int i6 = i5 + 1; i6 < nt; i6++)
							x_expand_3_6[n3_6++].Add3(t[i4], t[i5], t[i6]);
					}
				}
			}
		}
	}

	goto next;
	// going back, for a non empty index, count it back
back:
	if (--s3 >= spb3)goto next;
	// and set last index
	xindex3[nxindex3++].Open(0, n3_6, tcells);

}
void GCHK::DebugExpand() {
	cout << "debugexpand nindex=" << nxindex3 << " ntot=" << n3_6 << endl;
	for (int i = 0; i < nxindex3 - 1; i++) {
		XINDEX3 w = xindex3[i];
		int ideb = w.ideb, iend = xindex3[i + 1].ideb;
		cout << Char27out(w.cellsbf) << "\t" << w.ideb << endl;
		for (int i = 0; i < 3; i++) {
			int c = w.tcells[i];
			if (!(w.cellsbf & (1 << c))) cout << "erreur i=" << i << " c=" << c << endl;
		}
	}
}
void GCHK::GoBanda() {// processing band a
	// loop on index 3 
	for (int i3 = 0; i3 < nxindex3; i3++) {
		wi3 = xindex3[i3];
		Init3clues();
		cout << Char27out(wi3.cellsbf) << " i3=" << i3 << endl;
		//if (i3 > 10)continue;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		uint32_t indf = xindex3[i3 + 1].ideb;
		for (uint32_t i6 = wi3.ideb; i6 < indf; i6++) {
			if (i6 > 27638) return;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			if (!(i6 == 16340 || i6 == 23874 || i6 == 27638)) continue;
			wi3_6 = x_expand_3_6[i6];
			if (Init3_6clues()) {
				cout << Char27out(filt32) << "\t i6=" << i6 << endl;
				p_cpt1g[maxcluesb - sbb.ncrit]++;
				p_cpt1g[ncluesa +7]++;
				sbb.Go();
			}
			if (gchk.aigstop)return;
		}
	}
}

void GCHK::Init3clues() {
	p_cpt2g[2]++;
	memcpy(tclues, wi3.tcells, sizeof wi3.tcells);
	uint64_t *t = genuasb12.tua;
	uint32_t n = genuasb12.nua;
	{ // reduce the UA table
		memset(ntuasmini, 0, sizeof ntuasmini);
		ntua = 0;
		register uint64_t F = wi3.cellsbf;
		for (uint32_t i = 0; i < n; i++) {
			register uint64_t U = t[i] & BIT_SET_2X, Ua, Ub;
			if (U&F)continue; // UA hit
			Ub = U >> 32; 
			Ua = U & BIT_SET_27; 
			uint32_t countb = _popcnt32((uint32_t)Ub);
			if (countb > 3) {// store in in ua/ub mode
				Ub <<= 32;
				if (ntua < 512)tua[ntua++] = Ua + Ub;
				continue;
			}
			// this is a mini row UA in band b use sub table
			for (int im = 0, mask = 7; im < 9; im++, mask <<= 3) {// find the mini row
				if (!(Ub&mask)) continue;
				if (countb == 3) im += 27;// now index for the subtable
				else {// must be count 2 and in mini row
					register uint32_t bit = (uint32_t)Ub ^ mask;
					bitscanforward(im, bit);
				}
				if (ntuasmini[im] < 100)
					tuasmini[im][ntuasmini[im]++] = (uint32_t)Ua;
			}
		}
	}
	nactivemini = 0;
	for (int i = 0; i < 36; i++)if (ntuasmini[i])
		activemini[nactivemini++] = i;
}
int GCHK::Init3_6clues() {
	moreb.Init();
	nmoreof = nmoreif = 0;// init more outfield table

	// build tclues band a (3 clues done)
	register uint32_t w = wi3_6.d;
	ncluesa = (w & 7) + 3;
	if (ncluesa > 5) return 0; //<<<<<<<<<<<<<<<<<<<<<<<<<<<<< provisoire 
	p_cpt2g[3]++;
	if (ncluesa > 3) {
		for (int i = 3; i < ncluesa; i++) {
			w >>= 8;
			tclues[i] = w & 0xff;
		}
	}
	{// compute min max clues band B exit if no room for 19
		mincluesb = ncluesa;
		if (!(mode_equal & 1))mincluesb++;
		if (mincluesb == 7) return 0; //no 19 here
		mincluesc = mincluesb;
		if (!(mode_equal & 2))mincluesb++;
		maxcluesc=19- mincluesb- ncluesa;
		maxcluesb= 19 - ncluesa;
		if (!(mode_equal & 2))maxcluesb--;
		maxcluesb >>= 1;//to share in B and C
	}

	filt32 = wi3.cellsbf;// build the filter filt32
	filt32 |= wi3_6.cellsbf;
	{ // find active minirows in band b
		register uint32_t F = filt32;
		sbb.pairs27 = sbb.critbf = sbb.triplet = 0;
		for (uint32_t i = 0; i < nactivemini; i++) {
			uint32_t iac = activemini[i], *t = tuasmini[iac], n = ntuasmini[iac];
			for (uint32_t iu = 0; iu < n; iu++) {
				if (!(t[iu] & F)) {// Ua not hit, must be in band b
					if (iac < 27) sbb.pairs27 |= 1 << iac;
					else sbb.triplet |= 1 << (iac - 27);
					break; //first unsolved is ok
				}
			}
		}
	}
	{ // set up mini rows status and field
		sbb.mini1 = sbb.mini2 = sbb.mini3 = 0;
		register uint32_t P = sbb.pairs27;
		for (int im = 0, mask = 7, bit = 1; im < 9; im++, mask <<= 3, bit <<= 1) {// check mini rows
			register uint32_t M = P & mask;
			if (!M) {
				if (sbb.triplet&bit)sbb.critbf |= mask;
			}
			else {
				sbb.triplet &= ~bit;// triplet would be redundant
				uint32_t cc = _popcnt32(M);
				if (cc > 1) {
					sbb.critbf |= mask;
					if (cc == 2)sbb.mini2 |= bit;
					else sbb.mini3 |= bit;
				}
				else {
					sbb.critbf |= mask ^ M;
					sbb.mini1 |= bit;

				}
			}
		}
	}
	sbb.Ncrit();
	if (sbb.ncrit > maxcluesb)return 0;
	sbb.tuaif = btuaif;
	sbb.tuaof = btuaof;
	sbb.andoutf = BIT_SET_27;
	{	// build remaining uas in and out field
		sbb.nuaif = sbb.nuaof = 0;
		register uint32_t F = filt32, IF = sbb.critbf;
		for (uint32_t i = 0; i < ntua; i++) {
			register uint64_t U = tua[i];
			if ((uint32_t)U & F) continue;
			U >>= 32;// now ua bandb
			register uint32_t Ub = (uint32_t)U;
			if (Ub & IF) 				sbb.AddIF(Ub);
			else {
				if (sbb.ncrit == maxcluesb)return 0;
				sbb.AddOF(Ub);
				sbb.andoutf &=Ub;
				if ((!sbb.andoutf) && sbb.ncrit == (maxcluesb - 1))
					return 0;
			}
		}
	}
	{	// add band b uas
		uint32_t *t = bands_abc[1].tua, n = bands_abc[1].nua;
		register uint32_t  IF = sbb.critbf;
		for (uint32_t i = 0; i < n; i++) {
			register uint32_t Ub = t[i];
			if (Ub & IF) 		sbb.AddIF(Ub);
			else {
				if (sbb.ncrit == maxcluesb)return 0;
				sbb.AddOF(Ub);
				sbb.andoutf &= Ub;
				if ((!sbb.andoutf) && sbb.ncrit == (maxcluesb - 1))return 0;
			}
		}
	}
	nbif = sbb.nuaif;
	nbof = sbb.nuaof;
	{	//_______________ guas filter => already killed/forced plus sub table
		// build new subtables still active not fixed
		register uint32_t F = filt32;
		nguared_2 = 0;
		memset(ntuar2, 0, sizeof ntuar2);
		forced81_2.SetAll_0(); forced81_3.SetAll_0();
		for (int i = 0; i < genb12.nactive2; i++) {
			int i81 = genb12.tactive2[i];
			GEN_BANDES_12::SGUA2 & w = genb12.tsgua2[i81];
			uint32_t * tuasw = tuar2[i81], nt = 0;
			for (uint32_t iua = 0; iua < w.nua; iua++) {
				register uint64_t U = w.tua[iua] & BIT_SET_2X;
				register uint32_t Ua, Ub;
				Ub = U >> 32; Ua = U & BIT_SET_27; 
				if (Ua&F) continue; // hit by banda
				if (!Ub) {// empty in band 2 stop and force
					forced81_2.Set_c(i81);
					nt = 0;
					goto nexti81_2;
				}
				Ub |= (_popcnt32(Ub) << 27);
				AddUA32(tuasw, nt, Ub);
			}
			if (nt) {
				ntuar2[i81] = nt;
				guar2i81[nguared_2++] = i81;
			}
		nexti81_2:;
		}
		nguared_3 = 0;
		memset(ntuar3, 0, sizeof ntuar3);
		for (int i = 0; i < genb12.nactive3; i++) {
			int i81 = genb12.tactive3[i];
			GEN_BANDES_12::SGUA3 & w = genb12.tsgua3[i81];
			uint32_t * tuasw = tuar3[i81], nt = 0;
			for (uint32_t iua = 0; iua < w.nua; iua++) {
				register uint64_t U = w.tua[iua] & BIT_SET_2X;
				register uint32_t Ua, Ub;
				 Ub = U >> 32; Ua = U & BIT_SET_27;
				if (Ua&F) continue; // hit by banda
				if (!Ub) {// empty in band 2 stop and force
					forced81_3.Set_c(i81);
					nt = 0;
					goto nexti81_3;
				}
				Ub |= (_popcnt32(Ub) << 27);
				AddUA32(tuasw, nt, Ub);
			}
			if (nt) {
				ntuar3[i81] = nt;
				guar3i81[nguared_3++] = i81;
			}
		nexti81_3:;
		}
	}
	return 1;// can continue searching bandb solutions
}

/* results as example (34620 3_6)
130			0 max outfield b 0
1186		1 max outfield b 1
5074		2 max outfield b 2
9276		3 max outfield b 3
10171		4 max outfield b 4
6653		5 max outfield b 5
1886		6 max outfield b 6
51			7 max outfield b 7

78		11 band A 4 clues
2924		12 band A 5 clues
31425		13 band A 6 clues

*/


void BANDB::Go() {//start band b expansion
	if (diagbug) {
		cout << "start B expansion nmiss=" << nmiss
			<< " nuaif=" << nuaif << " nuaof=" << nuaof << "\t p_cpt2g[3]=" << p_cpt2g[3] << endl;
	}
	diag = 0;
	known_bb = rknown_bb = 0;
	ndead = BIT_SET_27;
	active_bb = critbf;// active cells in field
	wactive0 = BIT_SET_27 ^ active_bb;//  active cells out field
	{// compute minband b and loop from minb to maxb
		if (gchk.ncluesa < 5) {
			cout << gchk.ncluesa << " " << gchk.mincluesb << " " << gchk.maxcluesb
				<< " a; minb; maxb nuaof="<<nuaof << "\tandoutf="<<andoutf<<"\tncrit="<<ncrit ;
		}
		else return;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		int minof = 0, minw=gchk.mincluesb;
		if (nuaof)minof = (andoutf) ? 1 : 2;
		int minw2 = ncrit + minof;
		if (minw2 > minw)minw = minw2;// should never be > maxcluesb
		cout << "\t minw  " << minw << endl;
		for (ncluesb = minw; ncluesb <= gchk.maxcluesb; ncluesb++) {
			if (ncluesb < 6) 				p_cpt1g[15]++;
			else if (ncluesb >7) p_cpt1g[18]++;
			else  p_cpt1g[10+ncluesb]++;
			if (ncluesb != 7) continue;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			nmiss= ncluesb - ncrit;
			if (gchk.nmoreof) {// insert here moreof if any
				for (uint32_t i = 0; i < gchk.nmoreof; i++) {
					andoutf &= gchk.more_of[i];
					AddOF(gchk.more_of[i]);
				}
			}
			Status();
			BANDB hn = *this;	
			if (!nmiss) {if(!nuaof)hn.Go_Critical();}
			else if (nmiss == 1){//nothing to do if andoutf=0
				if( andoutf) hn.Go_miss1();
			}
			else if (nmiss == 2) hn.Go_miss2();
			else hn.ExpandBandB();
		}
	}

	//else Go_Not_Critical_missn();
}


void BANDB::Status() {
	cout << "Bandb Status ncrit=" <<ncrit
		<< " nuaif"<<nuaif<<" nuaof"<<nuaof
		<< " ncluesb=" << ncluesb << " nmiss=" << nmiss<< endl;
	cout << Char27out(critbf) << " in field bf" << endl ;
	cout << Char27out(pairs27) << " pairs27" << endl ;
	cout << Char9out(mini_all) << " all minis2" << endl;
	cout << Char9out(mini1) << "     minis1" << endl;
	cout << Char9out(mini2) << "     minis2" << endl;
	cout << Char9out(mini3) << "     minis3" << endl ;
	cout << Char9out(triplet) << " mini triplets" << endl << endl;

}
void BANDB::DebugIfOf() {
	cout << "In Field Bandb table nuaif= "<<nuaif << endl;
	for (uint32_t i = 0; i < nuaif; i++)
		cout << Char27out(tuaif[i]) << endl;
	cout << "Out Field Bandb table nuaof="<< nuaof << endl;
	for (uint32_t i = 0; i < nuaof; i++)
		cout << Char27out(tuaof[i]) << endl;

}
//======================== start band b expansion

//_____________________ critical
void BANDB::CriticalAssignCell(int Ru) {// assign a cell within the critical cells
	// Ru is usually a register containing a 27 bits field with one bit on
	// 2 pairs in a miniriow have already been applied
	known_bb |= Ru;
	uint32_t cell;
	bitscanforward(cell, Ru); // catch the cell
	register int mini = C_minirow[cell],// minirow to clear
		bit = 1 << mini,
		Mask = 7 << (3 * mini);
	if (bit & mini3) {// cell is in minirow with 3 pairs active
		active_bb &= ~Ru; //clear the cell
		mini3 ^= bit; // now only a pair to hit
		mini1 |= bit;
	}
	else {// either one pair or a triplet in the minirow
		//active_b3 &= ~Ru; //clear the cell
		active_bb &= (~Mask); // kill the minirow as active
		mini1 &= ~bit;
		triplet &= ~bit;
	}
}
int BANDB::IsFinalOrMultiple(uint32_t * wua) {
	// all if uas must be hit
	register int bf = known_bb | active_bb;
	for (uint32_t iua = 0; iua < nuaif; iua++) {
		register int Ru = tuaif[iua];
		if(! (Ru & bf)) return 1;//known UA not hit multiple 
	}
	if (gchk.moreb.Check(bf)) return 1;// more hit
	if (!active_bb) {
		gchk.CriticalFinalCheck(known_bb,diag);
		return 1; 
	}
	if (bf != rknown_bb) {
		rknown_bb = bf;
		if (gchk.IsMultiple(rknown_bb, diag)) return 1;
	}
	return 0;
}
int BANDB::BuildIF_short() {
	tuaif = gchk.tuaif;
	nuaif = 0;
	for (uint32_t iua = 0; iua < gchk.nbif; iua++) {
		register int Ru = gchk.btuaif[iua];
		if (Ru & known_bb) continue;// already hit, forget it
		Ru &= active_bb;
		if (!Ru) return 1;// dead branch
		Ru |= _popcnt32(Ru) << 27;
		AddUA32(tuaif, nuaif, Ru);
	}
	// also more if pending
	for (uint32_t iua = 0; iua < gchk.nmoreif; iua++) {
		register int Ru = gchk.more_if[iua];
		if (Ru & known_bb) continue;// already hit, forget it
		Ru &= active_bb;
		if (!Ru) return 1;// dead branch
		Ru |= _popcnt32(Ru) << 27;
		AddUA32(tuaif, nuaif, Ru);
	}
	return 0;
}

int BANDB::ShrinkUas1() {
	irloop = 0;
	uint32_t * tn = &tuaif[nuaif], n = 0;
	for (uint32_t iua = 0; iua < nuaif; iua++) {
		register int Ru = tuaif[iua];
		if (Ru & known_bb) continue;// already hit, forget it
		Ru &= active_bb;
		if (!Ru) return 1;// dead branch
		if (_popcnt32(Ru) == 1) {// assign it and reduce the active cells
			CriticalAssignCell(Ru);
			irloop = 1;// should loop for new singles
		}
		else {
			Ru |= _popcnt32(Ru) << 27;
			AddUA32(tn, n, Ru);
		}
	}
	tuaif = tn;
	nuaif = n;
	if (!nuaif) irloop = 0;// no need to loop again
	return 0;
}
void BANDB::Go_Critical(uint32_t * wua) {// critical situation all clues in pairs tripl:ets
	if (++p_cpt1g[20]<10) {
		cout << Char27out(active_bb) << "B active bb entry critical" << endl;
		cout << Char32out(known_bb) << "known 32 at entry" << endl;
		//diag = 1;
	}
	//else return;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	if (mini2) {// assign common cell 
		for (int i = 0, bit = 1, mask = 7; i < 9; i++ , bit <<= 1, mask <<= 3) {
			if (mini2&bit) {
				active_bb &= (~mask); // clear the minirow
				known_bb |= mask & (~pairs27);// and set the common cell a			
			}
		}
	}
	mini2 = 0;
	if (BuildIF_short()) return;// shrink and sort the table in field
	if (ShrinkUas1()) return;//assign possibles
	if (IsFinalOrMultiple(wua)) return;
	if (irloop)		CriticalLoop();
	else CriticalExitLoop();
}
void BANDB::CriticalLoop() {
	if (ShrinkUas1()) return;
	if (irloop)CriticalLoop();
	else CriticalExitLoop();
}
void BANDB::CriticalExitLoop() {
	int nmissb = ncluesb - _popcnt32(known_bb);// missing clues
	if ( diag) {
		cout << Char27out(known_bb) << "B entry critical exit loop nmissb="<< nmissb << endl;
	}
	if (nmissb < 0)return;
	if (IsFinalOrMultiple())return;
	if (nuaif) {		// find the smallest ua and apply it
		int wua = 0, sizeua = 27;
		uint32_t cell;
		if (nmissb == 1) {//most frequent case
			register int and_uas = active_bb;
			for (uint32_t i = 0; i < nuaif; i++) {
				and_uas &= tuaif[i];
			}
			if (!and_uas) return; // no possibility
			wua = and_uas;
		}
		else if (mini1) {	// use in priority an unsolved pair it is a smallest
			uint32_t mini;
			bitscanforward(mini, mini1);
			int  shift = 3 * mini, mask = 7 << shift;
			wua = active_bb & mask;// catch the minirow
		}
		else {
			for (uint32_t i = 0; i < nuaif; i++) {
				register int ua = tuaif[i]& active_bb,
					cc = _popcnt32(ua);
				if (cc < sizeua) { wua = ua; sizeua = cc; }
				if (cc < 3)break; // this is the minimum
			}
			if (sizeua >= 2 && triplet) {// use the triplet in priority
				uint32_t mini;
				bitscanforward(mini, triplet);
				int  shift = 3 * mini, mask = 7 << shift;
				wua = active_bb & mask;// catch the minirow
			}
		}
		if ( diag)		cout << Char27out(wua) << " wua to use" << endl;

		while (bitscanforward(cell, wua)) {
			if (diag)		cout << Char27out(wua) << " wua entry while nmissb=" << nmissb<< endl;
			register int bit = 1 << cell;
			wua ^= bit;// clear bit			
			active_bb ^= bit;//  now a dead cell downstream
			BANDB hn = *this;
			hn.CriticalAssignCell(bit);
			hn.CriticalLoop();
		}
	}
	else Critical_0_UA(); // no more ua, some non assigned pairs or triplets
}
void BANDB::Critical_0_UA() {
	int nmissb = ncluesb - _popcnt32(known_bb);// missing clues
	if (nmissb < 0)return;
	if (!nmissb) {// nothing more to assign (granted at first call in a branch)
		gchk.CriticalFinalCheck(known_bb,diag);
		return;
	}
	if (mini3) {// in active minirows with 3 pairs, assign 2
		while (mini3) {
			uint32_t mini;
			bitscanforward(mini, mini3);
			int shift = 3 * mini, bit = 1 << shift;
			mini3 ^= 1 << mini; //clear bit the mini row is always killed
			active_bb &= ~(7 << shift); // clear also the bitfield of active cells
			int tp[3][2] = { { 0, 1 }, { 0, 2 }, { 1, 2 } };
			for (int i = 0; i < 3; i++) {
				int * tpi = tp[i];
				BANDB hn = *this;
				hn.CriticalAssignCell(bit << tpi[0]);
				hn.CriticalAssignCell(bit << tpi[1]);
				hn.Critical_0_UA();
			}
		}
		return;
	}
	if (mini1) {// active pair in minirow
		uint32_t mini;
		bitscanforward(mini, mini1);
		int  shift = 3 * mini, bit = 1 << shift, mask = 7 << shift;
		int x = active_bb & mask;// catch the minirow
		active_bb &= ~mask;// clear the minirow
		mini1 ^= 1 << mini;// and clear the minirow bit as active
		for (int i = 0; i < 3; i++) {
			int bb = bit << i;
			if (x&bb) {
				BANDB hn = *this;
				hn.CriticalAssignCell(bb);
				hn.Critical_0_UA();
			}
		}
		return;
	}
	// now must be active triplet in minirow
	if (triplet) {// safety control should always be
		uint32_t mini;
		bitscanforward(mini, triplet);
		int shift = 3 * mini, bit = 1 << shift, mask = 7 << shift;
		active_bb &= ~mask;// clear the minirow
		triplet ^= 1 << mini;// and clear the minirow bit as active
		for (int i = 0; i < 3; i++) {
			int bb = bit << i;
			BANDB hn = *this;
			hn.CriticalAssignCell(bb);
			hn.Critical_0_UA();
		}
	}
}
// ____________________ sub critical
void  BANDB::Go_Subcritical() {
	//if (++p_cpt1g[21] < 30){
	//	cout << Char27out(known_bb) << "B entry sub critical nmiss= " << nmiss << endl;
	//	diag = 1;
	//}
	//else return;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	active_bb = active_sub = critbf;
	// check first if a global solution  is still possible  
	if (IsFinalOrMultiple())		return;// not valid using all cells
	uint32_t cct = _popcnt32(critbf) - ncrit;
	if (cct < nmiss)return;// not enough remaining cells in GUA2s GUA3s to reach the count
	ndead = 0;
	Go_SubcriticalMiniRow();// find the first miss
}
void BANDB::Go_SubcriticalMiniRow() {
	int c2[3] = { 3, 5, 6 };// 2 cells in a mini row
	int bit = 1 << ndead, mask = 7 << (3 * ndead);
	for (int i = ndead; i < 9; i++, bit <<= 1, mask <<= 3) {
		register int M = active_sub & mask;
		if (!M)continue;
		ndead = i;
		if (bit & mini1) {// gua2 pair assign both
			BANDB hn = *this;
			hn.mini1 ^= bit;
			hn.SubMini(M, mask);
		}
		else if (bit & mini2) {// 2 gua2 pairs assign 2 out of 3
			for (int j = 0; j < 3; j++) {
				BANDB hn = *this;
				hn.mini2 ^= bit;
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				hn.SubMini(M, mask);
			}
		}
		else if (bit & mini3) {// 3 gua2 pairs assign all
			BANDB hn = *this;
			hn.mini3 ^= bit;
			hn.SubMini(M, mask);
		}
		else if (bit & triplet) {// gua3 assign 2 out of 3
			for (int j = 0; j < 3; j++) {
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				BANDB hn = *this;
				hn.triplet ^= bit;
				hn.SubMini(M, mask);
			}
		}
		else { // second add in the mini row one residual cell take it
			BANDB hn = *this;
			hn.SubMini(M, mask);
		}
	}
}
void BANDB::SubMini(int M, int mask) {
	known_bb |= M;// assign 1 or 2
	nmiss--;// one added
	active_bb &= ~mask;
	active_sub ^= M;
	if (nmiss) Go_SubcriticalMiniRow();// continue till no missing clue 
	else 		Go_Critical();	// leave sub critical   enter  critical
}
//_____________________ not critical
void BANDB::ShrinkUasOf() {
	if (known_bb) {// shrink the out field table
		uint32_t * tn = &tuaof[nuaof], n = 0;
		register uint32_t Ra = wactive0,
			Rfilt = known_bb;
		andoutf = BIT_SET_27;
		for (uint32_t iua = 0; iua < nuaof; iua++) {
			register int Ru = tuaof[iua];
			if (Ru & Rfilt) continue;
			Ru &= Ra;
			if (!Ru)return; // empty need at least one outfield
			andoutf &= Ru;
			Ru |= _popcnt32(Ru) << 27;
			AddUA32(tn, n, Ru);
		}
		for (uint32_t i = 0; i < gchk.nmoreof; i++) {// insert moreof 
			register int Ru = gchk.more_of[i];
			if (Ru & Rfilt) continue;
			Ru &= Ra;
			if (!Ru)return; // empty need at least one outfield
			andoutf &= Ru;
			Ru |= _popcnt32(Ru) << 27;
			AddUA32(tn, n, Ru);
		}
		tuaof = tn;
		nuaof = n;
	}
}

void BANDB::Go_miss1() {// not called if more than 1 needed
	ShrinkUasOf();// if known from up stream
	if (1) {
		cout << Char27out(known_bb )<< " nmiss1 nuaof=" << nuaof << endl;
	}
	if (!nuaof) {// subcritical in hn
		BANDB hn = *this;
		hn.Go_Subcritical();// to test later
	}
	uint32_t wua = wactive0;
	if (nuaof) wua &= andoutf;
	if (1) {
		cout << Char27out(wua) << " nmiss1 wua"  << endl;
	}
	if (wua) { // apply first UA to use or all out field cells
		uint32_t res;
		while (bitscanforward(res, wua)) {
			int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
			BANDB hn = *this;	hn.nmiss--;	hn.known_bb |= bit;
			if (ncrit)hn.Go_Critical(&wua);
			else {// this is a final band b to test
				int uabr = gchk.IsMultiple(hn.known_bb, g17b.diag);
				if (uabr) {// multiple try to apply it upstream
					wua &= gchk.myuab;
				}
				else gchk.CriticalFinalCheck(hn.known_bb);
			}
		}
	}
}
void BANDB::Go_miss2() {
	ShrinkUasOf();// if known from up stream
	if (!nuaof) {// subcritical in hn
		BANDB hn = *this;
		hn.Go_Subcritical();// to test later		
	}
	uint32_t wua = wactive0;
	if (nuaof) {
		if (andoutf) wua &= andoutf;// one common to all uas
		else wua &= tua[0];// use first ua  
	}

	if (wua) { // apply first UA to use or all out field cells
		uint32_t res;
		while (bitscanforward(res, wua)) {
			int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
			BANDB hn = *this;	hn.nmiss--;	hn.known_bb |= bit;
			hn.Go_miss1();
		}
	}
}

void BANDB::ExpandBandB() {
	uint32_t nex = nmiss;
	struct SPB3 {// spots to find band 3 minimum valid solutions
		uint32_t  possible_cells, all_previous_cells, active_cells, iuabb;
	}spb3[10], *s3, *sn3;
	uint64_t ispot;
	s3 = spb3;
	s3->all_previous_cells = 0;
	s3->active_cells = wactive0;// all cells of active
	s3->iuabb = 0; // copy the start table
	s3->possible_cells = tuaof[0] & wactive0;
	int tcells[15];
	//____________________  here start the search
next:
	ispot = s3 - spb3;
	uint32_t cell;
	if (!bitscanforward(cell, s3->possible_cells))goto back;
	{// apply cell in bitfields
		register int bit = 1 << cell;
		tcells[ispot] = cell;
		s3->possible_cells ^= bit;// clear bit
		register int filter = s3->all_previous_cells | bit,
			ac = s3->active_cells ^ bit;
		sn3 = s3 + 1;
		sn3->all_previous_cells = filter;
		sn3->active_cells = s3->active_cells = ac;
		// nextspot:take the next available ua to loop
		for (uint32_t i = s3->iuabb + 1; i < nuaof; i++) {
			if (tuaof[i] & filter)continue;
			if (ispot >= nmiss - 1)goto next;
			sn3->iuabb = i;
			sn3->possible_cells = tua[i] & ac;
			s3 = sn3; // switch to next spot

			goto next;
		}
	}	// no more ua
	nmiss = ncluesb - (uint32_t)ispot - 1 - ncrit;
	if (ncrit) {
		if (!nmiss) {
			BANDB hn = *this;	hn.known_bb = sn3->all_previous_cells;
			hn.Go_Critical();
			goto next;
		}
	}
	// no in field subcritical or  add one outfield
	BANDB hn = *this;	hn.known_bb = sn3->all_previous_cells;
	hn.Go_missx();
	goto next;
back:
	if (--s3 >= spb3)goto next;
}
void BANDB::Go_missx() {// after no more ua in bandb expand
	{// subcritical in hn
		BANDB hn = *this;
		hn.Go_Subcritical(); 		
	}
	uint32_t wua = wactive0;
	if (wua) { // apply first UA to use or all out field cells
		uint32_t res;
		while (bitscanforward(res, wua)) {
			int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
			BANDB hn = *this;	hn.nmiss--;	hn.known_bb |= bit;
			if (hn.nmiss) hn.Go_missx();
			else hn.Go_Critical();
		}
	}
}

void BANDB::Go_Not_Critical_missn() {
	ShrinkUasOf();// if known from up stream
	if (!nuaof) {// subcritical in hn
		BANDB hn = *this;
		//hn.Go_Subcritical();// to test later		
	}
	uint32_t wua = wactive0;
	if (nuaof) {
		if (andoutf) wua &= andoutf;// one common to all uas
		else wua &= tua[0];// use first ua  
	}
	if (wua) { // apply first UA to use or all out field cells
		uint32_t res;
		while (bitscanforward(res, wua)) {
			int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
			BANDB hn = *this;	hn.nmiss--;	hn.known_bb |= bit;
			if (hn.nmiss > 2) hn.Go_Not_Critical_missn();
			else hn.Go_miss2();
		}
	}
}

//_________________ 
int GCHK::IsMultiple(int bf,int diag) {
	if (_popcnt32(bf) > 25) return 0;	
	nclues = ncluesa;// buil bandb part of tclues
	uint32_t cellb, wbf=bf;
	while( wbf){
		bitscanforward(cellb, wbf);
		wbf ^= 1 << cellb;
		tclues[nclues++] = cellb + 27;
	}
	p_cpt2g[7]++;
	register uint64_t myua = zh2b[0].ValidXY(tclues,nclues, 0);
	if (myua) {//store the fresh ua bands ab
		if (diag) cout << Char2Xout(myua) << "uaret band b multiple" << endl;
		register uint64_t uab = myua>>32,uaa= myua&BIT_SET_27,uab12=myua;
		uint64_t cc64 = _popcnt64(uab12);
		if (cc64 < 12) {// this should never be check for a bug 
			cout <<Char2Xout(uab12)<<"ua < 12 to add ??? g17b.npuz"<< g17b.npuz 
				<< endl;
			//sbb.Status();
			//sbb.DebugIfOf();
			//zh2b[0].DebugValidXY(tclues, nclues, 0);
		}
		if (cc64 < 20)genuasb12.AddUACheck(uab12 | ((uint64_t)cc64 << 59));// and update the UA table
		myuab = (uint32_t)uab;
		int cc = _popcnt32((uint32_t)uab);
		// add uab in field or in "more outfield
		uint32_t ua = myuab | _popcnt32(myuab) << 27;
		if (uab & sbb.critbf) {
			if (nmoreif < 128)AddUA32(more_if, nmoreif, ua);
		}
		else if (nmoreof < 128)AddUA32(more_of, nmoreof, ua);

		gchk.moreb.Add((uint32_t)uab);// see final or multiple

		if (cc < 2)return 1;// should never be <2
		if (cc >= 4) {// add it to the reduced table in ua ub mode 
			if (ntua < 512)		tua[ntua++] = myua;
			return 1;
		}
		int i36;
		for(int imini=0,bit=1,mask=7;imini<9;imini++,bit<<=1,mask<<=3){
			if(!(uab&mask) )continue;
			if (cc == 2) {// fresh mini2
				uint32_t bit27 = (uab&mask) ^ mask;
				bitscanforward(i36, bit27);
			}
			else i36 = imini+27;
			if (!ntuasmini[i36])activemini[nactivemini++] = i36;
			if (ntuasmini[i36] < 50)
				tuasmini[i36][ntuasmini[i36]++] =(uint32_t) uaa;
		}
	
	}
	return (myua>0);
}
void GCHK::GuasCollect(int bf) {

	final81_2 = forced81_2;// collect guas2 active
	for (uint32_t i = 0; i < nguared_2; i++) {
		int i81 = guar2i81[i];
		uint32_t * tua = tuar2[i81];
		for (uint32_t iua = 0; iua < ntuar2[i81]; iua++) {
			register uint32_t Ru = tua[iua];
			if (Ru&bf)continue;
			final81_2.Set_c(i81);
			break;// one ua not hit is enough here
		}
	}
	final81_3 = forced81_3;// collect guas3 active
	for (uint32_t i = 0; i < nguared_3; i++) {
		int i81 = guar3i81[i];
		uint32_t * tua = tuar3[i81];
		for (uint32_t iua = 0; iua < ntuar3[i81]; iua++) {
			register uint32_t Ru = tua[iua];
			if (Ru&bf)continue;
			final81_3.Set_c(i81);
			break;// one ua not hit is enough here
		}
	}
}
int GCHK::SetUpGuas2_3() {
	maxcluesc = 19 - nclues;
	GUAs & myb = gchk.guas;
	BF128 ws2 = final81_2 & myb.isguasocket2;
	BF128 ws3 = final81_3 & myb.isguasocket3;
	// switch to mini rows patterns
	int tix[81], ntix = ws2.Table3X27(tix);
	mini_bf1 = mini_bf2 = mini_bf3 = pairsbf = pairs27 = mini_triplet = 0;
	for (int i = 0; i < ntix; i++) {// switch from i_81 to mini rows
		int i81 = tix[i], imini = myb.ua2_imini[i81], bit = 1 << imini;
		if (mini_bf2&bit) mini_bf3 |= bit;
		if (mini_bf1&bit) mini_bf2 |= bit;
		mini_bf1 |= bit;
		pairsbf |= myb.ua_pair[i81];
		pairs27 |= (1 << myb.ua2_i27[i81]);
	}
	ntix = ws3.Table3X27(tix);// now triplets to mini rows
	for (int i = 0; i < ntix; i++) {// switch from i_81 to mini rows
		int imini = myb.ua3_imini[tix[i]], bit = 1 << imini;
		mini_triplet |= bit;
	}
	//___________________________ prepare a new band to process
	all_used_minis = mini_bf1 | mini_triplet;
	mini_triplet &= ~mini_bf1;// count only triplets with no pair
	mincount = _popcnt32(all_used_minis) + _popcnt32(mini_bf3);
	if (mincount > maxcluesc) return 1;// too many clues
	nmiss = ncluesb3 - mincount;
	mini_bf1 &= ~mini_bf2;// now pure one pair
	mini_bf2 &= ~mini_bf3;// now pure 2 pairs
	// set up pair + triplet bitfield
	if (mini_triplet) {// must add triplet minirow
		for (int i = 0, bit = 1, field = 7; i < 9; i++, bit <<= 1, field <<= 3)
			if (mini_triplet&bit)
				pairsbf |= field;
	}
	return 0;
}

void GCHK::EndCollectBand3() {
	//============= collect Gua46 and uas b3 for the band split them "in-field" "out-field"
	nuasb3_1 = nuasb3_2 = 0;
	andb3 = BIT_SET_27;
	register int  Rfilt = pairsbf;
	// first GUA46 usually shorter than UAs band3
	BF128  socket4 =guas.isguasocket4;// i81 3X
	socket4 &= final81_2;
	int * ua_46 = guas.ua_pair; // ua pattern
	int i81;
	while ((i81 = socket4.getFirstCell()) >= 0) {
		socket4.Clear_c(i81);// clear bit
		register uint32_t Ru = ua_46[i81] & BIT_SET_27;
		Ru |= _popcnt32(Ru) << 27;
		if (Ru & Rfilt)	AddUA32(uasb3_1, nuasb3_1, Ru);
		else {
			andb3 &= Ru;
			AddUA32(uasb3_2, nuasb3_2, Ru);
		}
	}
	uint32_t * to = bands_abc[2].tua;
	for (uint32_t i = 0; i < bands_abc[2].nua; i++) {
		register uint32_t Ru = to[i] & BIT_SET_27;
		Ru |= _popcnt32(Ru) << 27;
		if (Ru & Rfilt)	AddUA32(uasb3_1, nuasb3_1, Ru);
		else {
			andb3 &= Ru;
			AddUA32(uasb3_2, nuasb3_2, Ru);
		}
	}
}
void GCHK::CriticalFinalCheck(int bf,  int debug) {// no more ua is it a valid solution
	p_cpt2g[8]++;
	GuasCollect(bf);
	if (SetUpGuas2_3())return;
	if (debug)		cout << Char27out(bf) << "critical check bandb" << endl;
	register int ir = IsMultiple(bf,1);
	if (!ir) {// one valid bands A + B
		if (debug) 	cout << "band 2 valid for band 3" << endl;
		p_cpt2g[9]++;
		if (zhou[0].PartialInitSearch17(tclues, nclues))
			return;// would be  bug
		GoBand3();
	}
}
void GCHK::GoBand3() {
	if (gchk.aigstop)return;
	if (1) {
		char ws[82];
		strcpy(ws, empty_puzzle);
		ws[54] = 0;
		for (int i = 0; i < nclues; i++) {
			int cell = tclues[i];
			ws[cell] = zsol[cell];
		}
		fout1 << ws << "bab" << endl;
		return;
	}
	EndCollectBand3();
	mincluesc = maxcluesc - 2;
	uint32_t mini = mincount,mini2=nclues-ncluesa;
	if (!(mode_equal & 2))mini2++;// must be clues c >= clues b 
	if (nuasb3_2) {// mincount + 1 or more
		mini++;
		if (!andb3)mini++;
	}
	if (mini2 > mini)mini = mini2;
	if (mini > mincluesc)mincluesc = mini;
	for (ncluesb3 = mincluesc; ncluesb3 <= maxcluesc; ncluesb3++) {
		nmiss = ncluesb3 - mincount;
		if (nmiss < 0)continue;
		p_cpt2g[15] ++;
		//if (p_cpt2g[9] < 5)return;
		G17B3HANDLER hh0; hh0.Init();
		hh0.diagh = sbb.diag;
		if (!nmiss)hh0.Go_Critical();
		else if (nmiss == 1) {//nothing to do if andoutf=0
			if (andb3) hh0.Go_miss1b3();
		}
		else if (nmiss == 2) hh0.Go_miss2b3();
		//if (!mincount) ExpandBand3();
		//else 
		//else hh0.Go_Not_Critical_missn();
	}
}
int  GCHK::BuildUasB3_in(uint32_t known, uint32_t field) {
	nuas_in = 0;
	for (uint32_t iua = 0; iua < nuasb3_1; iua++) {
		register uint32_t  Ru = uasb3_1[iua];
		if (Ru & known) continue;// already hit, forget it
		Ru &= field;
		if (!Ru) return 1;// dead branch
		Ru |= _popcnt32(Ru) << 27;
		AddUA32(uas_in, nuas_in, Ru);
	}
	return 0;
}

void GCHK::Status() {
	cout << "Band3 Status"  << endl;
	cout << Char27out(pairsbf) << " pairs bf" << endl;
	cout << Char27out(pairs27) << " pairs 27" << endl;
	cout << Char9out(mini_bf1) << "     minis bf1" << endl;
	cout << Char9out(mini_bf2) << "     minis bf2" << endl;
	cout << Char9out(mini_bf3) << "     minis bf3" << endl;
	cout << Char9out(mini_triplet) << " mini triplets" << endl << endl;
}

//================ part 2  band 3 processing

int ZHOU::CallMultipleB3(ZHOU & o, uint32_t bf, int diag) {
	*this = o;
	if (diag) cout << Char27out(bf) << " call multipleb3" << endl;
	BF128 dca[9];
	int digitsbf = zh_g2.digitsbf;
	memcpy(dca, zh_g2.Digit_cell_Assigned, sizeof dca);
	{	
		uint32_t cc;
		register int x = bf;
		while (bitscanforward(cc, x)) {
			x ^= 1 << cc; //clear bit
			int cell = cc + 54, digit = gchk.grid0[cell];
			digitsbf |= 1 << digit;
			int xcell = cc + 64; // the cell value in 3x32 of a 128 bits map
			if (FD[digit][0].Off(xcell))  return 0;// check not valid entry
			Assign(digit, cell, xcell);
			dca[digit].Set(xcell);
		}
	}
	if (_popcnt32(digitsbf < 8)) {
		if (diag) {
			cout << "not eight digits" << endl;
			ImageCandidats();
		}
		return 1;// can not be one solution
	}
	zh_g2.s17_b3_mini = 1;
	BF128 w = cells_unsolved;
	w.bf.u32[3] = ~0;// keep rowunsolved settled
	for (int i = 0; i < 9; i++)  FD[i][0] &= w | dca[i];
	if (diag) ImageCandidats();
	//__________end assign last lot start solver
	zh_g.go_back = 0;	zh_g.nsol = 0; zh_g.lim = 1;// modevalid is set to  1
	int ir = Full17Update();
	if (diag) {
		cout << "after update" << endl;
		ImageCandidats();
	}
	if (ir == 2) return 0;// solved can not be multiple
	Guess17(0,diag);

	return zh_g.nsol;  
}
int ZHOU::Apply17SingleOrEmptyCells() {
	zh_g.single_applied = 0;
	// here  singles and empty cells till 4 cells searched 
	BF128 R1 = FD[0][0], R2 = R1 & FD[1][0]; 	R1 |= FD[1][0];
	BF128 R3 = R2 & FD[2][0]; R2 |= R1 & FD[2][0]; R1 |= FD[2][0];
	BF128 R4= R3 & FD[3][0]; 
		R3 |= R2 & FD[3][0]; R2 |= R1 & FD[3][0]; R1 |= FD[3][0];
	BF128 R5 = R4 & FD[4][0]; R4 |= R3 & FD[4][0];
		R3 |= R2 & FD[4][0]; R2 |= R1 & FD[4][0]; R1 |= FD[4][0];
	R5 |= R4 & FD[5][0]; R4 |= R3 & FD[5][0];
		R3 |= R2 & FD[5][0]; R2 |= R1 & FD[5][0]; R1 |= FD[5][0];
	R5 |= R4 & FD[5][6]; R4 |= R3 & FD[6][0];
		R3 |= R2 & FD[6][0]; R2 |= R1 & FD[6][0]; R1 |= FD[6][0];
	R5 |= R4 & FD[7][0]; R4 |= R3 & FD[7][0];
		R3 |= R2 & FD[7][0]; R2 |= R1 & FD[7][0]; R1 |= FD[7][0];
	R5 |= R4 & FD[8][0]; R4 |= R3 & FD[8][0];
		R3 |= R2 & FD[8][0]; R2 |= R1 & FD[8][0]; R1 |= FD[8][0];
	if ((cells_unsolved - R1).isNotEmpty()) 	return 1; // empty cells
	R1 -= R2; // now true singles
	R1 &= cells_unsolved; // these are new singles
	if (R1.isEmpty()) {// no single store pairs and more
		zh_g.pairs = R2 - R3;
		zh_g2.triplets = R3 - R4;
		zh_g2.quads = R4 - R5;
		return 0;
	}
	int tcells[80], ntcells = R1.Table3X27(tcells);
	for (int i = 0; i < ntcells; i++) {
		int cell = tcells[i];
		for (int idig = 0; idig < 9; idig++) {
			if (FD[idig][0].On_c(cell)) {
				Assign(idig, cell, C_To128[cell]);
				goto nextr1;
			}
		}
		return 1; // conflict with previous assign within this lot
	nextr1:;
	}
	zh_g.single_applied = 1;
	return 0;
}
int ZHOU::Full17Update() {
	if (zh_g.go_back) return 0;
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (Apply17SingleOrEmptyCells())	return 0; // locked empty cell or conflict singles in cells
		if (!zh_g.single_applied)	break;
	}
	return 1;
}
void ZHOU::Guess17(int index, int diag) {
	if (diag) {
		char ws[82];
		cout << zh_g.pairs.String3X(ws)<< " pairs "<< endl;
		cout << zh_g2.triplets.String3X(ws) << " triplets " << endl;
		cout << zh_g2.quads.String3X(ws) << " quads " << endl;
	}
	if (zh_g2.s17_b3_mini) {// look once for mini rows
		zh_g2.s17_b3_mini = 0;
		uint32_t b3 = cells_unsolved.bf.u32[2], aig = 0;
		if (!b3)return; // if b3 solved, all is solved
		for (uint32_t i = 0, mask = 7; i < 9; i++, mask <<= 3) {
			uint32_t mini = b3 & mask;
			if (_popcnt32(mini) < 2) continue;
			aig = 1;
			//cout << Char27out(mini) << " mini " << endl;
			// try this mini row as unsolved
			ZHOU * mynext = (this + 1);
			*mynext = *this;
			for (uint32_t j = 0, bit = 1; j < 27; j++, bit <<= 1) {
				if (bit&mini)continue;
				if (b3&bit) {
					int cell = j + 54, digit = zh_g2.grid0[cell];
					mynext->Seta_c(digit, cell);
				}
			}
			//cout << Char27out(mynext->cells_unsolved.bf.u32[2]) << " unsolved b3 " << endl;
			//cout << "appel suite zh_g.go_back="<< zh_g.go_back << endl;
			int ir = mynext->Full17Update();// solve as much as possible
			if (ir == 2) continue;
			uint32_t b3_n = mynext->cells_unsolved.bf.u32[2];
			if (!b3_n) continue; // now solved
			mynext->Guess17(0, 0);
			if (zh_g.nsol) return;
			zh_g.go_back = 0;// see why it is 1
		}
		if (aig) {
			int ir = Apply17SingleOrEmptyCells();// restore zh_g
			zh_g.go_back = 0;// see why it is 1
		}
	}
	BF128 w = zh_g.pairs;
	if (w.isEmpty())w = zh_g2.triplets;
	if (w.isEmpty())w = zh_g2.quads;
	if (w.isEmpty())w = cells_unsolved;
	{ // select band with more unsolved cells 
		uint32_t nfreecells = 0, nw;
		if (w.bf.u32[0]) {
			nfreecells = _popcnt32(cells_unsolved.bf.u32[0]);
		}
		if (w.bf.u32[1]) {
			if (nfreecells) {
				nw = _popcnt32(cells_unsolved.bf.u32[1]);
				if (nw > nfreecells) {
					nfreecells = nw;
					w.bf.u32[0] = 0;
				}
			}
			else	nfreecells = _popcnt32(cells_unsolved.bf.u32[1]);
		}
		if (w.bf.u32[2]) {
			if (nfreecells) {
				nw = _popcnt32(cells_unsolved.bf.u32[2]);
				if (nw > nfreecells) {
					nfreecells = nw;
					w.bf.u32[0] = 0;
					w.bf.u32[1] = 0;
				}
			}
		}
	}
	int xcell = w.getFirst128(),
		cell = From_128_To_81[xcell],
		digit = zh_g2.grid0[cell],
		tdig[10], ndig = 0;
	if (diag) {
		cout << "guess17 index=" << index << " cell " << cellsFixedData[cell].pt << endl;
	}
	// if first step try first false
	if(!index)	ClearCandidate_c(digit, cell);// force false
	for (int idig = 0; idig < 9; idig++)
		if (FD[idig][0].On(xcell))tdig[ndig++] = idig;
	for (int idig = 0; idig < ndig; idig++) {
		ZHOU * mynext = (this + 1);
		*mynext = *this;
		mynext->SetaCom(tdig[idig], cell, xcell);
		if (diag)cout <<"guess index="<<index<<" "
			<< tdig[idig]+1<< cellsFixedData[cell].pt << endl;
		mynext->Compute17Next(index+1,diag);
		if (zh_g.go_back) return;

	}
	if (!index) {
		FD[digit]->Set_c(cell);// restore the candidate
		SetaCom(digit, cell, xcell);
		if (diag)cout << "guess last index=" << index << " "
			<< digit + 1 << cellsFixedData[cell].pt << endl;
		Compute17Next(index, diag);

	}
}

void ZHOU::Compute17Next(int index, int diag) {
	int ir = Full17Update();
	if (!ir) return;// locked 
	if (diag>1) { cout << "index=" << index << endl; ImageCandidats(); }
	if (ir == 2) {//solved
		if (index) {// store false as ua
			BF128 & wua = zh_g2.cells_assigned;
			int * sol = gchk.grid0;
			wua.SetAll_0();;
			for (int i = 0; i < 81; i++) {
				int d = sol[i];
				if (FD[d][0].Off_c(i))	wua.Set_c(i);
			}
			zh_g.nsol++;
		}
		zh_g.go_back = 1;// closed anyway
		return;
	}
	Guess17(index , diag);// continue the process
}

void G17B3HANDLER::Init() {
	GCHK  & bab = gchk;
	nb3 = bab.ncluesb3;
	uasb3of = bab.uasb3_2;
	nuasb3of = bab.nuasb3_2;
	andoutf = bab.andb3;
	uasb3if = bab.uasb3_1;
	nuasb3if = bab.nuasb3_1;
	int tbitsrows[8] = { 0, 07, 07000, 07007, 07000000, 07000007, 07007000, 07007007 };
	known_b3 = rknown_b3 = 0;
	ndead = BIT_SET_27;
	active_b3 =pairsbf= bab.pairsbf;// active cells in field
	pairs27 = bab.pairs27;
	wactive0 = BIT_SET_27 ^ active_b3;//  active cells out field
	nmiss = bab.nmiss;
	ncritical = bab.mincount;
	nb3 = nmiss + ncritical;
	mini_bf1 = bab.mini_bf1;
	mini_bf2 = bab.mini_bf2;
	mini_bf3 = bab.mini_bf3;
	mini_triplet = bab.mini_triplet;
	diagh = 0;
}

uint32_t G17B3HANDLER::IsMultiple(int bf) {
	if (bf == rknown_b3) return 0;
	if (_popcnt32(bf) > 25) return 0;
	uint32_t ua = 0;
	rknown_b3 = bf;
	GCHK & bab = gchk;
	// check first if all tuab3 is hit
	int ir = zhou[1].CallMultipleB3(zhou[0], bf, 0);
	if (ir) {//consider store the fresh ua b3
		BF128 wua = zh_g2.cells_assigned;
		ua = wua.bf.u32[2];
		int cc = _popcnt32(ua);
		if ( diagh) cout << Char27out(ua) << " ua b3 retour" << endl;
		//if (1) cout << Char27out(ua) << " ua b3 retour" << endl;
		if ( cc < 4) {// needs more tests on performance evolution
			if (cc == 2) {// fresh gua2
				int i81 = gchk.GetI81_2(ua);
				if (i81 >= 0) {
					bab.final81_2.Set_c(i81);// new valid gua2 for other bands
					if (!bab.ntuar2[i81]) bab.guar2i81[bab.nguared_2++] = i81;
					if (bab.ntuar2[i81] < GUAREDSIZE)
						bab.tuar2[i81][bab.ntuar2[i81]++] = wua.bf.u32[1];
					GEN_BANDES_12::SGUA2 & sg = genb12.tsgua2[i81];
					if (sg.nua < SIZETGUA) 
						AddUA64(sg.tua, sg.nua, wua.bf.u64[0]);
				}
			}
			if (cc == 3) {// fresh gua2
				int i81 = gchk.GetI81_3(ua);
				if (i81 >= 0) {
					bab.final81_3.Set_c(i81);// new valid gua2 for other bands
					if (!bab.ntuar3[i81]) bab.guar3i81[bab.nguared_3++] = i81;
					if (bab.ntuar3[i81] < GUAREDSIZE)
						bab.tuar3[i81][bab.ntuar3[i81]++] = wua.bf.u32[1];
					GEN_BANDES_12::SGUA3 & sg = genb12.tsgua3[i81];
					if (sg.nua < SIZETGUA)
						AddUA64(sg.tua, sg.nua, wua.bf.u64[0]);
				}
			}
		}
	}
	return ua;
}

void G17B3HANDLER::Critical2pairs() {// assign 2 pairs in minirow to common cell
	int tbitsrows[8] = { 0, 07, 07000, 07007, 07000000, 07000007, 07007000, 07007007 };
	if (mini_bf2) {// and 2 pairs in minirow forced to common cell
		register int Rst = 07007007;// stack 0 pattern 
		for (int ist = 0; ist < 3; ist++) {
			int shrink = TblShrinkMask[mini_bf2 & (0111 << ist)];
			if (shrink) {// minirows 2 pairs in that stack
				register int Mask = tbitsrows[shrink] << (3 * ist);
				active_b3 &= (~Mask); // clear the minirow
				known_b3 |= Mask & (~pairs27);// and set the common cell as assigned
			}
		}
		mini_bf2 = 0;
	}
}
//================= critical process
void G17B3HANDLER::CriticalAssignCell(int Ru){// assign a cell within the critical cells
	// Ru is usually a regidster containing a 27 bits field with one bit on
	// 2 pairs in a miniriow have already been applied
	known_b3 |= Ru;
	uint32_t cell;
	bitscanforward(cell, Ru); // catch the cell
	register int mini = C_minirow[cell],// minirow to clear
		bit = 1 << mini,
		Mask = 7 << (3 * mini);
	if (bit & mini_bf3){// the cell is in a minirow with 3 pairs active
		active_b3 &= ~Ru; //clear the cell
		mini_bf3 ^= bit; // now only a pairto hit
		mini_bf1 |= bit;
	}
	else{// either one pair or a triplet in the minirow
		active_b3 &= (~Mask); // kill the minirow as active
		mini_bf1 &= ~bit;
		mini_triplet &= ~bit;
	}
}
int G17B3HANDLER::BuildIfShortB3() {
	if (gchk.BuildUasB3_in(known_b3, active_b3))return 1;
	uasb3if = gchk.uas_in;
	nuasb3if = gchk.nuas_in;
	return 0;
}

int G17B3HANDLER::ShrinkUas1() {
	irloop = 0;
	uint32_t * tn = &uasb3if[nuasb3if], n = 0;
	register uint32_t Ra = active_b3,
		Rfilt = known_b3;
	for (uint32_t iua = 0; iua < nuasb3if; iua++) {
		register int Ru = uasb3if[iua];
		if (Ru & known_b3) continue;// already hit, forget it
		Ru &= active_b3;
		if (!Ru) return 1;// dead branch
		if (_popcnt32(Ru) == 1) {// assign it and reduce the active cells
			CriticalAssignCell(Ru);
			irloop = 1;// should loop for new singles
		}
		else {
			Ru |= _popcnt32(Ru) << 27;
			AddUA32(tn, n, Ru);
		}
	}
	uasb3if = tn;
	nuasb3if = n;
	if (!n ) irloop = 0;// no need to loop again
	return 0;

}
void G17B3HANDLER::Go_Critical(){// critical situation all clues in pairs tripl:ets
	//if (g17b.debug17 > 1 && known_b3)cout << Char27out(known_b3) << " entry critical" << endl;
	if (1) {
		cout << Char27out(known_b3) << "entry critical " << endl;
		return;
	}
	p_cpt2g[16]++;
	active_b3 = pairsbf;
	Critical2pairs();// assign 2 pairs in minirow to common cell
	if (!active_b3){
		CriticalFinalCheck();
		return; // should be seen earlier if possible
	}
	if (IsMultiple(known_b3 | active_b3)) return;
	
	if(BuildIfShortB3())return;
	if (ShrinkUas1()) return;// dead branch
	if (!active_b3) {
		CriticalFinalCheck();
		return; // should be seen earlier if possible
	}
	if (IsMultiple(known_b3 | active_b3))	return;// not valid using all cells
	if (irloop)		CriticalLoop();
	else CriticalExitLoop();
}

void G17B3HANDLER::CriticalLoop(){
	if (ShrinkUas1()) return;
	if (irloop)CriticalLoop();
	else CriticalExitLoop();
}
void G17B3HANDLER::CriticalExitLoop(){
	int nmissb = gchk.ncluesb3 - _popcnt32(known_b3);// missing clues
	if (nmissb < 0)return;
	if (!active_b3){// nothing more to assign 
		if (nmissb)return;// dead cell in a mini row 3 pairs
		CriticalFinalCheck();
		return;
	}
	// check known + active with brute force
	int wknown = known_b3 | active_b3;
	if (IsMultiple(known_b3 | active_b3))	return;// not valid using all cells
	if (g17b.diag > 2|| diagh) cout << "critical exit loop moreto do nuasb3=" << nuasb3if << endl;

	if (nuasb3if){		// find the smallest ua and apply it
		int wua = 0, sizeua = 27;
		uint32_t cell;
		if (nmissb == 1) {//most frequent case
			register int and_uas = active_b3;
			for (uint32_t i = 0; i < nuasb3if; i++) {
				and_uas &= uasb3if[i];
			}
			if (!and_uas) return; // no possibility
			wua = and_uas;
		}
		else if (mini_bf1) {	// use in priority an unsolved pair it is a smallest
			uint32_t mini;
			bitscanforward(mini, mini_bf1);
			int  shift = 3 * mini, mask = 7 << shift;
			wua = active_b3 & mask;// catch the minirow
		}
		else {
			for (uint32_t i = 0; i < nuasb3if; i++) {
				register int ua = uasb3if[i]& active_b3,
					cc = _popcnt32(ua);
				if (cc < sizeua) { wua = ua; sizeua = cc; }
				if (cc < 3)break; // this is the minimum
			}
			if (sizeua >= 2 && mini_triplet) {// use the triplet in priority
				uint32_t mini;
				bitscanforward(mini, mini_triplet);
				int  shift = 3 * mini, mask = 7 << shift;
				wua = active_b3 & mask;// catch the minirow

			}
		}

		while (bitscanforward(cell, wua)){
			register int bit = 1 << cell;
			wua ^= bit;// clear bit
			// clean the bit in active_b3, this is now a dead cell downstream
			active_b3 ^= bit;
			G17B3HANDLER hn = *this;
			hn.CriticalAssignCell(bit);
			hn.CriticalLoop();
		}
	}
	else Critical_0_UA(); // no more ua, some non assigned pairs or triplets
}
void G17B3HANDLER::Critical_0_UA(){
	int nmissb = gchk.ncluesb3 - _popcnt32(known_b3);// missing clues
	if (g17b.diag >= 2) {
		cout << Char27out(known_b3) << "critical_0ua missb=" << nmissb << endl;
		PrintStatus();
	}
	if (nmissb < 0)return;
	if (!nmissb) {// nothing more to assign (granted at first call in a branch)
		CriticalFinalCheck();
		return;
	}
	if (mini_bf3) {// in active minirows with 3 pairs, assign 2
		while (mini_bf3) {
			uint32_t mini;
			bitscanforward(mini, mini_bf3);
			int shift = 3 * mini, bit = 1 << shift;
			mini_bf3 ^= 1 << mini; //clear bit the mini row is always killed
			active_b3 &= ~(7 << shift); // clear also the bitfield of active cells
			int tp[3][2] = { { 0, 1 }, { 0, 2 }, { 1, 2 } };
			for (int i = 0; i < 3; i++) {
				int * tpi = tp[i];
				G17B3HANDLER hn = *this;
				hn.CriticalAssignCell(bit << tpi[0]);
				hn.CriticalAssignCell(bit << tpi[1]);
				hn.Critical_0_UA();
			}
		}
		return;
	}
	if (mini_bf1) {// active pair in minirow
		uint32_t mini;
		bitscanforward(mini, mini_bf1);
		int  shift = 3 * mini, bit = 1 << shift, mask = 7 << shift;
		int x = active_b3 & mask;// catch the minirow
		active_b3 &= ~mask;// clear the minirow
		mini_bf1 ^= 1 << mini;// and clear the minirow bit as active
		for (int i = 0; i < 3; i++) {
			int bb = bit << i;
			if (x&bb) {
				G17B3HANDLER hn = *this;
				hn.CriticalAssignCell(bb);
				hn.Critical_0_UA();
			}
		}
		return;
	}
	// now must be active triplet in minirow
	if (mini_triplet) {// safety control should always be
		uint32_t mini;
		bitscanforward(mini, mini_triplet);
		int shift = 3 * mini, bit = 1 << shift, mask = 7 << shift;
		active_b3 &= ~mask;// clear the minirow
		mini_triplet ^= 1 << mini;// and clear the minirow bit as active
		for (int i = 0; i < 3; i++) {
			int bb = bit << i;
			G17B3HANDLER hn = *this;
			hn.CriticalAssignCell(bb);
			hn.Critical_0_UA();
		}
	}
}
void G17B3HANDLER::PrintStatus() {
	cout << "G17B3HANDLER Band3 Status" << endl;
	cout << Char27out(pairsbf) << " pairs bf" << endl;
	cout << Char27out(pairs27) << " pairs 27" << endl;
	cout << Char9out(mini_bf1) << "     minis bf1" << endl;
	cout << Char9out(mini_bf2) << "     minis bf2" << endl;
	cout << Char9out(mini_bf3) << "     minis bf3" << endl;
	cout << Char9out(mini_triplet) << " mini triplets" << endl << endl;
}
void G17B3HANDLER::CriticalFinalCheck(){// no more ua is it a valid solution 
	//if (p_cpt2g[17] == 2252)
	int ncl = _popcnt32(known_b3);
	//if (g17b.debug17) cout << "final check test ncl=" << ncl  << endl;
	if (ncl != gchk.ncluesb3) return; // should be seen earlier if possible
	register int ir = IsMultiple(known_b3);// , g17b.debug17);
	if (ir)		return;
	cout << "one sol to print final check valid id"<< p_cpt2g[15]
		<<Char32out(known_b3) << endl;
	char ws[82];
	strcpy(ws, empty_puzzle);
	for (int i = 0; i < gchk.nclues; i++) {
		int cell = gchk.tclues[i];
		ws[cell] = gchk.grid0[cell] + '1';
	}
	for (int i = 0, bit = 1; i < 27; i++, bit <<= 1)if (known_b3 & bit)
		ws[54 + i] = gchk.grid0[54 + i] + '1';
	fout1 << ws << ";"   << endl;
	p_cpt2g[25]++;
}
//=============== sub critical process   missing(s)  in the critical area
void G17B3HANDLER::Go_SubcriticalMiniRow() {
	int c2[3] = { 3, 5, 6 };// 2 cells in a mini row
	int bit = 1 << ndead, mask = 7 << (3 * ndead);
	for (int i = ndead; i < 9; i++, bit <<= 1, mask <<= 3) {
		register int M = active_sub & mask;
		if (!M)continue;
		ndead = i;
		if (g17b.diag >= 2)	cout << Char27out(M) << "  subcritical i= "<<i << endl;
		if (bit & mini_bf1) {// gua2 pair assign both
			G17B3HANDLER hn = *this;
			hn.mini_bf1 ^= bit;
			hn.SubMini(M, mask);
		}
		else if (bit & mini_bf2) {// 2 gua2 pairs assign 2 out of 3
			for (int j = 0; j < 3; j++) {
				G17B3HANDLER hn = *this;
				hn.mini_bf2 ^= bit;
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				hn.SubMini(M, mask);
			}
		}
		else if (bit & mini_bf3) {// 3 gua2 pairs assign all
			G17B3HANDLER hn = *this;
			hn.mini_bf3 ^= bit;
			hn.SubMini(M, mask);
		}
		else if (bit & mini_triplet) {// gua3 assign 2 out of 3
			for (int j = 0; j < 3; j++) {
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				G17B3HANDLER hn = *this;
				hn.mini_triplet ^= bit;
				hn.SubMini(M, mask);
			}
		}
		else { // second add in the mini row one residual cell take it
			G17B3HANDLER hn = *this;
			hn.SubMini(M, mask);
		}
	}
}
void G17B3HANDLER::SubMini( int M, int mask){
	if (g17b.diag > 2) {
		cout << Char27out(M)    << "entry subcritical M  nmiss="<<nmiss << endl;
		cout << Char27out(mask) << "                  mask " << endl;
		cout << Char27out(active_sub) << "                  active sub " << endl;
	}

	known_b3 |= M;// assign 1 or 2
	nmiss--;// one added
	active_b3 &= ~mask;
	active_sub ^= M;
	if (nmiss) Go_SubcriticalMiniRow();// continue till a"no missing clue condition"
	else {	// leave sub critical mode and enter the critical mode
		if (g17b.diag > 2) cout << "exit submini" << endl;
		Critical2pairs();// assign 2 pairs in minirow to common cell
//		if (g17b.bands_ab.BuildUas_in(known_b3, active_b3))return;
		if (ShrinkUas1() )return;// dead branch
		if (IsMultiple(known_b3 | active_b3))return;// not valid using all cells
		if (g17b.diag > 2)cout << "sub irloop=" << irloop<< " nuasb3=" << nuasb3if << endl;
		if (irloop)		CriticalLoop();
		else CriticalExitLoop();
	}
}
void G17B3HANDLER::Go_Subcritical(){// nmiss to select in the critical field
	if (g17b.diag >= 2) {
		cout << Char27out(known_b3) << "entry subcritical " << endl;
		cout << Char27out(active_b3) << "active b3 " << endl;
	}
	p_cpt2g[16]++;
	active_b3 = active_sub =pairsbf;
	// check first if a global solution  is still possible
	if (IsMultiple(known_b3 | active_b3))return;// not valid using all cells
	int cct = _popcnt32(pairsbf) - ncritical;
	if (g17b.diag > 2)	cout  << "subcritical cct="<<cct << endl;
	if (cct < nmiss)return;// not enough remaining cells in GUA2s GUA3s to reach the count
	ndead = 0;
	if (BuildIfShortB3())return;
	if (g17b.diag > 2)	cout << Char27out(active_sub) << "active sub subcritical " << endl;
	Go_SubcriticalMiniRow();// find the first miss
}
//======================================================================= not critical sequence
void G17B3HANDLER::ShrinkUasOfB3() {
	if (known_b3) {// shrink the out field table
		uint32_t * tn = &uasb3of[nuasb3of], n = 0;
		register uint32_t Ra = wactive0,
			Rfilt = known_b3;
		andoutf = BIT_SET_27;
		for (uint32_t iua = 0; iua < nuasb3of; iua++) {
			register int Ru = uasb3of[iua];
			if (Ru & Rfilt) continue;
			Ru &= Ra;
			if (!Ru)return; // empty need at least one outfield
			andoutf &= Ru;
			Ru |= _popcnt32(Ru) << 27;
			AddUA32(tn, n, Ru);
		}
		uasb3of = tn;
		nuasb3of = n;
	}
}

void G17B3HANDLER::Go_miss1b3() {// not called if more than 1 needed
	ShrinkUasOfB3();// if known from up stream
	if (1) {
		cout << Char27out(known_b3) << " nmiss1b3 nuaof=" << nuasb3of << endl;
	}
	if (!nuasb3of) {// subcritical in hn
		G17B3HANDLER hn = *this;
		hn.Go_Subcritical();// to test later
	}
	uint32_t wua = wactive0;
	if (nuasb3of) wua &= andoutf;
	if (1) {
		cout << Char27out(wua) << " nmiss1b3 wua" << endl;
	}
	if (wua) { // apply first UA to use or all out field cells
		uint32_t res;
		while (bitscanforward(res, wua)) {
			int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
			G17B3HANDLER hn = *this;	hn.nmiss--;	hn.known_b3 |= bit;
			if (ncritical)hn.Go_Critical();
			else {// this is a final band b to test
				int uabr = gchk.IsMultiple(hn.known_b3, diagh);
				if (uabr) {// multiple try to apply it upstream
					wua &= gchk.myuab;
				}
				else gchk.CriticalFinalCheck(hn.known_b3);
			}
		}
	}
}
void G17B3HANDLER::Go_miss2b3() {
	ShrinkUasOfB3();// if known from up stream
	if (!nuasb3of) {// subcritical in hn
		G17B3HANDLER hn = *this;
		hn.Go_Subcritical();// to test later		
	}
	uint32_t wua = wactive0;
	if (nuasb3of) {
		if (andoutf) wua &= andoutf;// one common to all uas
		else wua &= uasb3of[0];// use first ua  
	}

	if (wua) { // apply first UA to use or all out field cells
		uint32_t res;
		while (bitscanforward(res, wua)) {
			int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
			G17B3HANDLER hn = *this;	hn.nmiss--;	hn.known_b3 |= bit;
			hn.Go_miss1b3();
		}
	}
}


void G17B3HANDLER::Go_Not_Critical_missn() {
	if ( diagh) 	cout<<Char27out(known_b3) << "entry not_critical miss " << nmiss
		<<" nuasb3of="<< nuasb3of << endl;
	if (known_b3) {// shrink the out field table
		uint32_t * tn = &uasb3of[nuasb3of], n = 0;
		register uint32_t Ra = wactive0,
			Rfilt = known_b3;
		andoutf = BIT_SET_27;
		for (uint32_t iua = 0; iua < nuasb3of; iua++) {
			register int Ru = uasb3of[iua];
			if (Ru & Rfilt) continue;
			Ru &= Ra;
			if (!Ru)return; // empty need at least one outfield
			andoutf &= Ru;
			Ru |= _popcnt32(Ru) << 27;
			AddUA32(tn, n, Ru);
		}
		uasb3of = tn;
		nuasb3of = n;
	}
	else if (nmiss == 1 ) {
		andoutf = wactive0;
		for (uint32_t iua = 0; iua < nuasb3of; iua++)
			andoutf &= uasb3of[iua];
	}
	uint32_t wua = andoutf;
	if (nmiss > 1) wua = uasb3of[0] & BIT_SET_27;
	if (!nuasb3of)wua = wactive0;
	if (g17b.debug17 > 1 || diagh)		cout << Char27out(wua) << "wua to use  " << endl;

	if(wua){ // apply first UA to use or all out field cells 
		uint32_t res;
		while (bitscanforward(res, wua)) {
			if (g17b.debug17 > 1|| diagh)		cout << Char27out(wua) << "wua used nmiss= "<<nmiss << endl;
			int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
			G17B3HANDLER hn = *this; hn.nmiss--; hn.known_b3 |= bit;
			if (hn.nmiss) {
				hn.Go_Not_Critical_missn();
			}
			else hn.Go_Critical();
			
		}
	}
	if (!nuasb3of)	Go_Subcritical();// finish in Subcritical if no ua
}

void GCHK::ExpandBand3() {
	struct SPB3 {// spots to find band 3 minimum valid solutions
		GINT64 stack_count;
		uint32_t  possible_cells, all_previous_cells, active_cells, iuab3;
	}spb3[10], *s3, *sn3;
	uint32_t * tua = uasb3_2, nua = nuasb3_2;
	uint64_t ispot;
	s3 = spb3;
	s3->all_previous_cells = 0;
	s3->active_cells = BIT_SET_27;// all cells active
	s3->iuab3 = 0; // copy the start table
	s3->possible_cells = tua[0] & BIT_SET_27;
	int tcells[15];
	//____________________  here start the search
next:
	ispot = s3 - spb3;
	uint32_t cell;
	if (!bitscanforward(cell, s3->possible_cells))goto back;
	{// apply cell in bitfields
		register int bit = 1 << cell;
		tcells[ispot] = cell;
		s3->possible_cells ^= bit;// clear bit
		register int filter = s3->all_previous_cells | bit,
			ac = s3->active_cells ^ bit;
		sn3 = s3 + 1;
		sn3->all_previous_cells = filter;
		sn3->active_cells = s3->active_cells = ac;
		// nextspot:take the next available ua to loop
		for (uint32_t i = s3->iuab3 + 1; i < nua; i++) {
			if (tua[i] & filter)continue;
			if (ispot >= ncluesb3 - 1)goto next;
			sn3->iuab3 = i;
			sn3->possible_cells = tua[i] & ac;
			s3 = sn3; // switch to next spot
			goto next;
		}
	}	// no more ua
	{	// check if this is a valid band 1+2+3 (can not be a valid 16)
		int ir = zhou[1].CallMultipleB3(zhou[0], sn3->all_previous_cells, 0);
		if (ir) {//consider store the fresh ua b3
			uint32_t ua = zh_g2.cells_assigned.bf.u32[2];
			if (nua < 500) // 500 is the limit for tuasb2
				tua[nua++] = ua;
			if (ispot < ncluesb3 - 1) {// if not a 17 do next
				sn3->iuab3 = nua - 1;
				sn3->possible_cells = ua & sn3->active_cells;
				s3 = sn3; // switch to next spot
			}
			goto next;
		}
		if (ispot < ncluesb3 - 1) {//  not a 17 should never be
			cout << " bug false 17" << endl;
			cerr << " bug false 17" << endl;
		}
		// valid 17
		cout << "one sol to print expand b3 valid id" << p_cpt2g[15] << endl;
		char ws[82];
		strcpy(ws, empty_puzzle);
		for (int i = 0; i < nclues; i++) {
			int cell = tclues[i];
			ws[cell] = gchk.grid0[cell] + '1';
		}
		for (int i = 0; i <= ispot; i++)
			ws[54 + tcells[i]] = gchk.grid0[54 + tcells[i]] + '1';
		fout1 << ws <<  ";expand" << endl;
		p_cpt2g[26]++;
		goto next;
	}
back:
	if (--s3 >= spb3)goto next;
}


//=============================== debugging sequences
void GCHK::GodebugInit(int mode) {
	cout << zsol << " traité" << endl;
	cout << "ua bands1+2   \t" << genuasb12.nua << endl;
	cout << "guas socket2  \t" << genb12.ntua2 << endl;
	cout << "guas socket3  \t" << genb12.ntua3 << endl;
	cout << "active socket2\t" << genb12.nactive2 << endl;
	cout << "active socket3\t" << genb12.nactive3 << endl;
	if (mode & 1) {
		cout << "table uas" << endl;
		uint64_t *t = genuasb12.tua;
		uint32_t n = genuasb12.nua;
		for (uint32_t i = 0; i < n; i++) cout << Char2Xout(t[i]) << endl;

	}
	if (mode & 2) {
		cout << "sockets 2 table" << endl;
		int n2 = 0;
		for (int i = 0; i < genb12.nactive2; i++) {
			int i81 = genb12.tactive2[i];
			GEN_BANDES_12::SGUA2 & w = genb12.tsgua2[i81];
			cout << i81 << " " << w.nua << endl;
			n2 += w.nua;
		}
		cout << "cumul=" << n2 << endl;
		cout << "sockets 3 table" << endl;
		int n3 = 0;
		for (int i = 0; i < genb12.nactive3; i++) {
			int i81 = genb12.tactive3[i];
			GEN_BANDES_12::SGUA3 & w = genb12.tsgua3[i81];
			cout << i81 << " " << w.nua << endl;
			n3 += w.nua;
		}
		cout << "cumul=" << n3 << " total=" << n2 + n3 << endl;
	}

}

