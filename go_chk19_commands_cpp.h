
//========================================
const char * libs_c17_00_cpt1g[40] = {
	"0 max outfield b 0",//0
	"1 max outfield b 1",//1
	"2 max outfield b 2",//2
	"3 max outfield b 3",//3
	"4 max outfield b 4",//4
	"5 max outfield b 5",//5
	"6 max outfield b 6",//6
	"7 max outfield b 7",//7
	"8 max outfield b 8",//8
	"9 ",//9
	"10 band A 3 clues",//10
	"11 band A 4 clues",//11
	"12 band A 5 clues",//12
	"13 band A 6 clues",//13
	"14 ",//14
	"15 band B <6 clues",//15
	"16 band B 6 clues",//16
	"17 band B 7 clues",//17
	"18 band B >7 clues",//18
	"19 ",//19
	"20 ",//20
	"21 ",//21
	"22 ",//22
	"23 ",//23
	"24 ",	"25 ",	"26 ",	"27",	"28", "29",
	"30",	"31",	"32",	"33", "34",
	"35",	"36",	"37",	"38", "39",

};
const char * libs_c17_00_cpt2g[40] = {
	"0 bands1+2 processed entries M10",//0
	"1 total bands 3",//1
	"2 index 3 entries",//2
	"3 index 3_5 entries",//3
	"4 bb crit ",//4
	"5 bb miss12 ",//5
	"6 bb miss more ",//6
	"7 calls brute force",//7
	"8 XY brute force",//8
	"9 valid brute force",//9
	"10 critical band3",//10
	"11 not critical 1",//11
	"12 not critical 234",//12
	"13 not critical 56",//13
	"14 ",//14
	"15 entry band3 handler excluding critical+ua outfield",//15
	"16 critical + sub critical",//16
	"17 add 1 from active",//17
	"18 n uas at start",//18
	"19 n gua2s at start  ",//19
	"20 n gua3s at start  ",//20
	"21 n sockets2",//21
	"22 n sockets1",//22
	"23 max bands 3",//23
	"24 max bands 3 go",//24
	"25 number of 17 found std",
	"26 number of 17 found through expandb3",
	"27",	"28", "29",
	"30",	"31",	"32",	"33", "34",
	"35",	"36",	"37",	"38", "39",

};
void Go_c17_00() {// p2 process
/*
	cout << "Go_c17_00 search batch 17 clues 656 566 " << endl;
	cout << sgo.vx[0] << " -v0- band 0_415" << endl;
	cout << sgo.vx[2] << " -v2- skip first nnn restart after batch failure" << endl;
	cout << sgo.vx[3] << " -v3- last entry number for this batch must be > vx[2]" << endl;
	cout << sgo.vx[4] << " -v4- 0 if p2a 1 if p2b" << endl;

	int it16_start = sgo.vx[0];
	g17b.debug17 = g17b.aigstop=0;
	g17b.diag = sgo.vx[6];
	genb12.skip = sgo.vx[2];
	genb12.last = sgo.vx[3];
	if (sgo.vx[2] < 0) {
		cerr << "invalid value for skip" << endl;
		return;
	}
	if (sgo.vx[3] < sgo.vx[2]) {
		cerr << "invalid value for last to process" << endl;
		return;
	}
	zh_g.modevalid = 1;
	zh_g2.grid0 = genb12.grid0;
	zh_g2.zsol = zh_g2.stdfirstsol;
	memset(p_cptg, 0, sizeof p_cptg);// used in debugging sequences only
	memset(p_cpt1g, 0, sizeof p_cpt1g);// used in debugging sequences only
	memset(p_cpt2g, 0, sizeof p_cpt2g);// used in debugging sequences only
	genb12.Start(0);
	genb12.NewBand1(sgo.vx[0]);
	cout << "print final stats" << endl;
	for (int i = 0; i < 40; i++) {
		if (!p_cpt2g[i])continue;
		cout << p_cpt2g[i] << "\t\t" << libs_c17_00_cpt2g[i] << endl;
	}
	cout << "exit final stats" << endl;
*/
}

//=========================entry file of solution grids to search
void Go_chk19_10() {
	STD_B416 bax[3];
	zh_g.modevalid = 1;
	zh_g2.grid0 = gchk.grid0;
	zh_g2.zsol = zh_g2.stdfirstsol;
	// search 17 using a file having known  as entry and one 17 given 6 6 5
	char * ze = finput.ze;
	int * zs0 = gchk.grid0, npuz = 0;
	cout << "Go_chk19_10() search 17-19 in a solution grid " << endl;
	while (finput.GetLigne()) {
		if(strlen(ze)<81) continue;// skip blank lines
		npuz++;
		memset(p_cpt1g, 0, sizeof p_cpt1g);
		gchk.npuz = npuz;
		gchk.aigstop= 0;
		long tdeb = GetTimeMillis();
		cout << ze <<  " to process  n="  << npuz << endl;

		// ====catch entry uas and rank
		for (int i = 0; i < 81; i++)zs0[i] = ze[i] - '1';
		BANDMINLEX::PERM perm_ret;
		bandminlex.Getmin(zs0, &perm_ret);
		int ib1 = perm_ret.i416, ib1a = t416_to_n6[ib1];
		bax[0].InitBand2_3(ib1, ze, perm_ret);
		cout << ib1 << "\tib1 mins=" << t16_min_clues[bax[0].i416] << endl;
		bandminlex.Getmin(&zs0[27], &perm_ret);
		int ib2 = perm_ret.i416, ib2a = t416_to_n6[ib2];
		bax[1].InitBand2_3(ib2, &ze[27], perm_ret);
		cout << ib2 << "\tib2 mins=" << t16_min_clues[bax[1].i416] << endl;
		bandminlex.Getmin(&zs0[54], &perm_ret);
		int ib3 = perm_ret.i416, ib3a = t416_to_n6[ib3];
		bax[2].InitBand2_3(ib3, &ze[54], perm_ret);
		cout << ib1  << "\t" << ib2  << "\t" << ib3  
			<<"\tmins="<< t16_min_clues[bax[0].i416] 
			<< t16_min_clues[bax[1].i416] 
			<< t16_min_clues[bax[2].i416] << endl;
		int tsort[3];
		{// sort entry increasing order of min clues
			for (int i = 0; i < 3; i++) 
				tsort[i] = i+(t16_min_clues[bax[i].i416]<<8);
			for(int i=0;i<2;i++) for(int j=i+1;j<3;j++)
				if (tsort[i] > tsort[j]) {
					int temp = tsort[i];
					tsort[i] = tsort[j];
					tsort[j]=temp;
				}
		}
		{// now start all perms in bands a b c mode 
			int tsort[3] = { 0,1,2 };// pour test uas sur 17
			int mode[6] = {3,2,1,0,1,0};// 
			for (int ip = 0; ip < 6; ip++) {// 6 bands perms
				int * tp = tperm6[ip];
				gchk.Start(bax, tp, tsort, mode[ip]);
				break;
			}
		}


		//genb12.nband3 = 1;
		//myband1.DoExpandBand();// expand band1
		//genb12.ValidInitGang();
		//g17b.npuz = npuz;
		//g17b.GoM10();
		cout << "print puzzle  stats" << endl;
		for (int i = 0; i < 20; i++) {
			if (!p_cpt1g[i])continue;
			cout << p_cpt1g[i] << "\t\t" << libs_c17_00_cpt1g[i] << endl;
		}

	}
	cout << "print final stats" << endl;
	for (int i = 0; i < 20; i++) {
		if (!p_cpt2g[i])continue;
		cout << p_cpt2g[i] << "\t\t" << libs_c17_00_cpt2g[i] << endl;
	}
}
/* processing abc increasing order tperm6
1 2=1 3=2
1 2=1 3>2  a<=b<=c
1 2>1 3=2
1 2>1 3>2

1 3=1 2>3  a<=b<c
1 3>1 2>3

2 1>2 3=1  a<b<=c
2 1>2 3>2 

2 3>2 1>3  a<b<c

3 1>3 2=1  a<b<=c
3 1>3 2>1

3 2>3 1>2 a<b<c


*/

