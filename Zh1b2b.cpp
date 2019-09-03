/*
Based on code posted to <http://forum.enjoysudoku.com/3-77us-solver-2-8g-cpu-testcase-17sodoku-t30470.html>
by user zhouyundong_2012.
The copyright is not specified.
*/
#include "main.h"   
#include "Zh1b2b.h"

// row unsolved 6 bits per digit lookup table
uint64_t zh2b_t_runsolved[9] = { 077 , 077 << 6 , 077 << 12 , 077 << 18 , 077 << 24 ,
(uint64_t)077 << 32 ,(uint64_t)077 << 38 ,(uint64_t)077 << 44 ,(uint64_t)077 << 50 };
uint32_t zh2b_t_runsolvedshift[9] = { 0,6,12,18,24,32,38,44,50 };
//========================= global variable and working areas 
extern ZH_GLOBAL zh_g;
ZH2B zh2b[40]; //  
//ZH2B zh2b_i, zh2b_i1;
ZH2B_GLOBAL   zh2b_g;   // 2 bands 9 digits
ZH2B5_GLOBAL   zh2b5_g;   // 2_5 digits 2 bands
ZH2B_1D_GLOBAL zh2b1d_g;  // one digit 2 bands

ZH2B_1D zh2b1d[6]; // maxi 6 guesses, one per row
ZH2B5 zh2b5[10]; // solved digit per digit

ZHONE_GLOBAL zh1b_g;
ZHONE zhone[20]; // one band 9 digits 
ZHONE zhone_i; 


//#include "go_17sol_zx_UaCollector_cpp.h"
//GENUAS_1Bx genuas1b;
//============================= ZH_GLOBAL code
ZH2B_GLOBAL::ZH2B_GLOBAL(){
	zsol = 0; // no solution unless required buy the user
	nctlg =  0;

}
void ZH2B_GLOBAL::GetBands(int * g1, int * g2) {
	// build sol per digit and pm per digit at start
	memset(fd_sols, 0, sizeof fd_sols);
	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 3; j++) {
			int cell = 9 * j + i, dig = puz0[cell];
			fd_sols[1][dig].bf.u32[0] |= Zhoucol << i;
			fd_sols[1][dig].bf.u32[1] |= Zhoucol << i;
			fd_sols[0][dig].bf.u32[0] |= 1 << cell;
			dig = puz0[cell+27];
			fd_sols[1][dig].bf.u32[0] |= Zhoucol << i;
			fd_sols[1][dig].bf.u32[1] |= Zhoucol << i;
			fd_sols[0][dig].bf.u32[1] |= 1 << cell;
		}
	}
}

void ZH2B_GLOBAL::InitGangster(int * g0, int * g1) {
	uint64_t col64 = (uint64_t)Zhoucol | ((uint64_t)Zhoucol << 32);
	memcpy(fd_revised, fd_sols[1], sizeof fd_revised);
	for (int i = 0; i < 9; i++, col64 <<= 1) {
		if (g0[i] == g1[i])continue;
		int changes = g0[i] ^ g1[i]; // one added one cleared
		for (int d = 0, bit = 1; d < 9; d++, bit <<= 1) {// check digits
			if (!(changes & bit)) continue;
			if (g0[i] & bit)fd_revised[d] &= ~col64;
			else fd_revised[d] |= col64;
		}		
	}
}


uint64_t ZH2B_GLOBAL::BuildUaret(BF64 * wsol) {
	ua_ret = 0;
	for (int i = 0; i < 9; i++) {
		BF64 w = wsol[i] - fd_sols[0][i];
		ua_ret |= w.bf.u64;
	}
	return ua_ret;
}


//============= zh2B code for uas 2 bands and validity 2 bands puzzle
/*
void ZH2b::Start_nFloors(int floors, BF64 * mypm) {
	cells_unsolved;bf.u64=BIT_SET_2X^ solved_cells;
	for (uint32_t i = 0; i < nd; i++) {
		myfd[i] &= cells_unsolved;
		if (myfd[i].count == 6)return ;
	}
	
}*/
//============================================ ZH2B code
void ZH2B::Init_std_bands() {//init after zh1b_g getband
	zh2b_g.ndigits = 9;
	memcpy(this, zh2b_start, sizeof zh2b_start);
	memcpy(FD, zh2b_g.fd_sols[1], sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
}
void ZH2B::Init_gang() {//init after zh1b_g InitGangster
	zh2b_g.ndigits = 9;
	memcpy(this, zh2b_start, sizeof zh2b_start);
	memcpy(FD, zh2b_g.fd_revised, sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
}

void ZH2B::DebugSol() {//init after zh1b_g getband
	zh2b_g.ndigits = 9;
	memset(this, 0, sizeof zh2b[0]);
	memcpy(FD, zh2b_g.fd_sols[0], sizeof FD);
	ImageCandidats();
}

#define UPDN(I,J)A=FD[I].bf.u32[J];\
Shrink = (TblShrinkMask[A & 0x1FF] | \
TblShrinkMask[ (A>>9) & 0x1FF]<<3 | \
TblShrinkMask[ (A>>18) & 0x1FF]<<6);\
if ((A &=TblComplexMask[Shrink]) ==0)  return 0; \
S = ((A | (A >> 9) | (A >> 18)) & 0x1FF); \
FD[I].bf.u32[1 - J] &= TblMaskSingle[S]; \
S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]]; \
CompFD[I].bf.u32[J] = FD[I].bf.u32[J] = A;

#define UPWCL(I,P,Q,R,T,U,V,W,X)cl = ~(A & TblRowMask[S]); \
cells_unsolved.bf.u32[I] &= cl; \
wcl[P]&= cl;wcl[Q]&= cl;wcl[R]&= cl;wcl[T]&= cl;\
wcl[U]&= cl;wcl[V]&= cl;wcl[W]&= cl;wcl[X]&= cl;

#define UPWCL5(I,P,Q,R,T)cl = ~(A & TblRowMask[S]); \
cells_unsolved.bf.u32[I] &= cl; \
wcl[P]&= cl;wcl[Q]&= cl;wcl[R]&= cl;wcl[T]&= cl;


#define UPDN1(I,P,Q,R,T,U,V,W,X)Shrink = (TblShrinkMask[A & 0x1FF] | \
TblShrinkMask[ (A>>9) & 0x1FF]<<3 | \
TblShrinkMask[ (A>>18) & 0x1FF]<<6);\
if ((A &=TblComplexMask[Shrink]) ==0)  return 0; \
S = ((A | (A >> 9) | (A >> 18)) & 0x1FF); \
S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]]; \
if ((A>>27) != S){\
cl = ~(A & TblRowMask[S]); \
A=(A&BIT_SET_27) | (S<<27);\
cells_unsolved &= cl; \
wcl[P]&= cl;wcl[Q]&= cl;wcl[R]&= cl;wcl[T]&= cl;\
wcl[U]&= cl;wcl[V]&= cl;wcl[W]&= cl;wcl[X]&= cl;}\
CompFD[I] = FD[I] = A



//================================ZH2B code

void ZH2B::InitTclues(uint32_t * tclues, int n) {
	memset(zh2b_g.Digit_cell_Assigned_init, 0, sizeof zh2b_g.Digit_cell_Assigned_init);
	for (int icell = 0; icell < n; icell++) {
		int cell = tclues[icell], digit = zh2b_g.puz0[cell];
		int xcell = C_To128[cell]; // the cell value in 3x32 of a 128 bits map

		Assign(digit, cell, xcell);
		zh2b_g.Digit_cell_Assigned_init[digit].Set(xcell);
	}
	for (int i = 0; i < 9; i++)  FD[i] &= cells_unsolved | 
		zh2b_g.Digit_cell_Assigned_init[i];
}

int ZH2B::ApplySingleOrEmptyCells() {// only singles
	zh2b_g.single_applied = 0;
	uint64_t * map = &FD[0].bf.u64;
	uint64_t unsolved = cells_unsolved.bf.u64;
	register uint64_t R2 = map[0] & map[1],
		R1 = (map[0] | map[1]), Map;// digits 12
	Map = map[2]; R2 |= R1 & Map; R1 |= Map;
	Map = map[3];  R2 |= R1 & Map; R1 |= Map;
	Map = map[4];  R2 |= R1 & Map; R1 |= Map;
	Map = map[5];  R2 |= R1 & Map; R1 |= Map;
	Map = map[6];  R2 |= R1 & Map; R1 |= Map;
	Map = map[7];  R2 |= R1 & Map; R1 |= Map;
	Map = map[8];  R2 |= R1 & Map; R1 |= Map;
	if (unsolved & (~R1)) return 1; // locked
	R1 &= ~R2;
	R1 &= unsolved; // forget solved seen as singles
	if (R1) zh2b_g.single_applied = 1;
	else return 0;
	while (R1) {// usually a very small number of cells to assign
		uint32_t res;
		if (!bitscanforward64(res, R1)) break;
		uint64_t bit = (uint64_t)1 << res; // switch to the bit value
		R1 &= ~bit;  // clear the bit
		// call Seta(int digit, int xcell) so find digit
		for (int idig = 0; idig < 9; idig++) {
			if (map[idig] & bit) {// this is the digit
				//				if (FD[idig].Off(res))  return 1; // invalid, gane locked
				//				Seta(idig, res);
				int cell = From_128_To_81[res];
				Assign(idig, cell, res);
				goto nextr1;// finished for that cell
			}
		}
		return 1; //conflict with a previous cell assugn
	nextr1: {}
	}
	return 0;// not locked
}
int ZH2B::Seta(int digit, int xcell) { // single in cell
	int cell = From_128_To_81[xcell],
		block = TblBoard_Block[cell];
	if (FD[digit].Off(xcell)) return 1; // not valid
	Assign(digit, cell, xcell);
	BF64 *Fd = &FD[digit];
	*Fd &= AssignMask_Digit[cell].u64[0];
	int ddig = 6 * digit;
	if (digit > 4) ddig += 2;// second bloc of 32 bits
	rows_unsolved.Clear(ddig + C_row[cell]);//6*digit + row
	cells_unsolved.Clear(xcell);
	BF64 * RF = &FD[8];
	for (; RF >= FD; RF--)RF->Clear(xcell);
	Fd->Set(xcell); // restore bit for digit assigned
	return 0;
}

int ZH2B::Update(){
	int Shrink = 1;
	register int S, A;
	register unsigned int cl, *wcl = FD[0].bf.u32;
	while (Shrink ){
	  Shrink = 0;
	  if (!rows_unsolved.bf.u32[0])goto digit5;

	  {register unsigned int  AR = rows_unsolved.bf.u32[0];// valid for digits 0,1,2,3,4
	  if (!(AR & 077))goto digit1;

//=digit 0
	  if (FD[0].bf.u32[0] ==  CompFD[0].bf.u32[0])goto digit0b;
		UPDN(0,0)	if ((AR & 7) != S){
				AR &= 07777777770 | S;	UPWCL(0, 2,4,6,8,10,12,14,16)	}

digit0b:if (FD[0].bf.u32[1] ==CompFD[0].bf.u32[1])goto digit1;
		UPDN(0,1)	if (((AR >> 3) & 7) != S){
				AR &= 07777777707 | (S << 3);	UPWCL(1, 3,5,7,9,11,13,15,17)	}

digit1:	if (!(AR  & 07700))goto digit2;

		if (FD[1].bf.u32[0] ==	CompFD[1].bf.u32[0])goto digit1b;
		UPDN(1,0)	if (((AR >> 6) & 7) != S){
			AR &= 07777777077 | (S << 6);UPWCL(0, 0, 4, 6, 8, 10, 12, 14, 16)	}

digit1b:if (FD[1].bf.u32[1] ==	CompFD[1].bf.u32[1])goto digit2;
		UPDN(1,1)		if (((AR >> 9) & 7) != S){
			AR &= 07777770777 | (S << 9);UPWCL(1, 1, 5, 7, 9, 11, 13, 15, 17)	}

digit2:	if (!(AR  & 0770000))goto digit3;

		if (FD[2].bf.u32[0] ==	CompFD[2].bf.u32[0])goto digit2b;
		UPDN(2,0)	if (((AR >> 12) & 7) != S){
			AR &= 07777707777 | (S << 12);	UPWCL(0, 0, 2, 6, 8, 10, 12, 14, 16)}

digit2b:if (FD[2].bf.u32[1] ==	CompFD[2].bf.u32[1])goto digit3;
		UPDN(2,1)	if (((AR >> 15) & 7) != S){
			AR &= 07777077777 | (S << 15);	UPWCL(1, 1, 3, 7, 9, 11, 13, 15, 17)	}

digit3: if (!(AR & 077000000))goto digit4;

		if (FD[3].bf.u32[0] == CompFD[3].bf.u32[0])goto digit3b;
		  UPDN(3,0)	  if (((AR >> 18) & 7) != S){
			  AR &= 07770777777 | (S << 18); UPWCL(0, 0, 2, 4, 8, 10, 12, 14, 16)	 }

digit3b:  if (FD[3].bf.u32[1] == CompFD[3].bf.u32[1])goto digit4;
		  UPDN(3, 1)if (((AR >> 21) & 7) != S){
			  AR &= 07707777777 | (S << 21); UPWCL(1, 1, 3, 5, 9, 11, 13, 15, 17)	}

digit4:if (!(AR & 07700000000))goto end01234;

		if (FD[4].bf.u32[0] ==	CompFD[4].bf.u32[0])goto digit4b;
		UPDN(4, 0)if (((AR >> 24) & 7) != S){
			AR &= 07077777777 | (S << 24);  UPWCL(0, 0, 2, 4, 6, 10, 12, 14, 16)	}

digit4b:if (FD[4].bf.u32[1] == CompFD[4].bf.u32[1])goto end01234;
		UPDN(4, 1)if (((AR >> 27) & 7) != S){
			AR &= 0777777777 | (S << 27);  UPWCL(1, 1, 3, 5, 7, 11, 13, 15, 17)	}

end01234: rows_unsolved.bf.u32[0] = AR;
		}// end of validity for AR 01234


digit5:
	  if (!rows_unsolved.bf.u32[1])continue; // second lot  4 digits

	  {register unsigned int  AR = rows_unsolved.bf.u32[1];// valid for digits 5,6,7,8
	  if (!(AR & 077))goto digit6;

	  if (FD[5].bf.u32[0] ==	CompFD[5].bf.u32[0])goto digit5b;
		UPDN(5, 0) if ((AR & 7) != S){
			AR &= 07777777770 | S;	UPWCL(0, 0, 2, 4, 6, 8, 12, 14, 16)		}

digit5b:if (FD[5].bf.u32[1] == CompFD[5].bf.u32[1])goto digit6;
		UPDN(5, 1) 	if (((AR >> 3) & 7) != S){
			AR &= 07777777707 | (S << 3); UPWCL(1, 1, 3, 5, 7, 9, 13, 15, 17)	}

digit6:  if (!(AR & 07700))goto digit7;

		if (FD[6].bf.u32[0] ==  CompFD[6].bf.u32[0])goto digit6b;
		UPDN(6, 0) if (((AR >> 6) & 7) != S){
			AR &= 07777777077 | (S << 6);	UPWCL(0, 0, 2, 4, 6, 8, 10, 14, 16)	  }

digit6b: if (FD[6].bf.u32[1] == CompFD[6].bf.u32[1])goto digit7;
		UPDN(6, 1) if (((AR >> 9) & 7) != S){
			AR &= 07777770777 | (S << 9); UPWCL(1, 1, 3, 5, 7, 9, 11, 15, 17)  }

digit7:   if (!(AR  & 0770000))goto digit8;

		if (FD[7].bf.u32[0] ==			  CompFD[7].bf.u32[0])goto digit7b;
		  UPDN(7, 0)if (((AR >> 12) & 7) != S){
			  AR &= 07777707777 | (S << 12);  UPWCL(0, 0, 2, 4, 6, 8, 10, 12, 16)  }

digit7b:  if (FD[7].bf.u32[1] ==  CompFD[7].bf.u32[1])goto digit8;
		  UPDN(7, 1) if (((AR >> 15) & 7) != S){
			  AR &= 07777077777 | (S << 15);  UPWCL(1, 1, 3, 5, 7, 9, 11, 13, 17)  }

digit8:  if (!(AR  & 077000000))goto end5678;

		  if (FD[8].bf.u32[0] ==  CompFD[8].bf.u32[0])goto digit8b;
		  UPDN(8,0)	  if (((AR >> 18) & 7) != S){
			  AR &= 07770777777 | (S << 18);  UPWCL(0, 0, 2, 4, 6, 8, 10, 12, 14)	  }

digit8b: if (FD[8].bf.u32[1] == CompFD[8].bf.u32[1])goto end5678;
		  UPDN(8,1)	  if (((AR >> 21) & 7) != S){
			  AR &= 07707777777 | (S << 21);  UPWCL(1, 1, 3, 5, 7, 9, 11, 13, 15)  }

end5678:rows_unsolved.bf.u32[1] = AR;
	  }// end of validity for AR
  }// end while
#ifdef DIAG
	cout << "end update cycle" << endl;
	Debug(1);
	ImageCandidats();
#endif
  return 1;
}
int ZH2B::FullUpdate() {
	if (zh2b_g.go_back) return 0;
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (ApplySingleOrEmptyCells())			return 0; // locked empty cell or conflict singles in cells
		if (zh2b_g.single_applied) {
			continue;
		}
		break;
	}
	return 1;
}

char * ZH2B::SetKnown(char * zs) {
	strcpy(zs, &empty_puzzle[27]);
	int tdig[9];// build the table of digits
	register int A = rows_unsolved.bf.u32[0];
	tdig[0] = A & 077; A >>= 6; tdig[1] = A & 077; A >>= 6; tdig[2] = A & 077; A >>= 6;
	tdig[3] = A & 077; A >>= 6; tdig[4] = A & 077;
	A = rows_unsolved.bf.u32[1];
	tdig[5] = A & 077; A >>= 6; tdig[6] = A & 077; A >>= 6;
	tdig[7] = A & 077; A >>= 6; tdig[8] = A & 077;

	for (int idig = 0; idig < 9; idig++) {// one digit
		int arowsj = tdig[idig];// 6 rows
		if (arowsj == 077) continue;
		for (int ib = 0; ib < 2; ib++) {// 3 blocs per digit
			int arows = (arowsj >> (3 * ib)) & 7;
			if (arows == 7) continue; // not assigned
			unsigned int band = FD[idig].bf.u32[ib];
			for (int j = 0; j < 3; j++) if (!(arows & (1 << j))) {
				int row = (band >> TblMult9[j]) & 0x1ff;
				uint32_t  irow;
				bitscanforward(irow, row);
				int	cell = Tblstartblock[ib] + TblMult9[j] + irow;
				zs[cell] = idig + '1';
			}
		}
	}
	return zs;
}
uint64_t ZH2B::ValidXY(uint32_t * tclues, int n, int test) {
	if (test)cout << "entry validxy" << endl;
	zh2b_g.ua_ret = 0;
	Init_std_bands();
	InitTclues(tclues, n);
	if (FullUpdate()) {
		if (test)	ImageCandidats();
		if (rows_unsolved.isEmpty()) return 0;// solved 
		if(n>15){ //try to find a UA minirow in band b
			register uint32_t b2 = rows_unsolved.bf.u32[1];
			if (!b2)return 0; // if b2 solved, all is solved
			for (uint32_t i = 0, mask = 7; i < 9; i++, mask <<= 3) {
				uint32_t mini = b2 & mask;
				if (_popcnt32(mini) < 2) continue;
				ZH2B * mynext = (this + 1);
				*mynext = *this;
				for (uint32_t j = 0, bit = 1; j < 27; j++, bit <<= 1) {
					if (bit&mini)continue;
					if (b2&bit) {
						int cell = j + 27, digit = zh2b_g.puz0[cell];
						mynext->SetaC(digit, cell);
					}
				}
				mynext->ComputeNext();
				if (zh2b_g.ua_ret) {
					if (test) {
						cout << Char2Xout(zh2b_g.ua_ret) << " ua mini " << i << endl;
					}

					return zh2b_g.ua_ret;// return if  ua found
				}
			}
		}
		// if does not work, suspect valid try best case first
		if (test) cout << "no mini row ua , try best" << endl;
		GuessValidB12();//std solve
	}
	return zh2b_g.ua_ret;
}
void ZH2B::DebugValidXY(uint32_t * tclues, int n, int test) {
	if (test)cout << "entry validxy" << endl;
	zh2b_g.ua_ret = 0;
	Init_std_bands();
	InitTclues(tclues, n);
	ImageCandidats();
	FullUpdate();
	ImageCandidats();
	for (int i = 0; i < 9; i++) 		FD[i] = zh2b_g.fd_sols[0][i];
	cout << "registered solution" << endl;
	ImageCandidats();
}
void ZH2B::GuessGo(int dig, BF64 & wsol) {// done in a new ocurrence
	*this = *(this - 1);	// apply wsol and make next step
	FD[dig] = wsol;// update will do the job
	ComputeNext();
}
void ZH2B::ComputeNext() {
	int ir = FullUpdate();
	if (ir == 1)GuessValidB12();
	else if (ir == 2) {// solved 
		zh2b_g.ua_ret = 0;
		for (int i = 0; i < 9; i++) {
			BF64 w = FD[i] - zh2b_g.fd_sols[0][i];
			zh2b_g.ua_ret |= w.bf.u64;
		}
	}
}
void ZH2B::GuessValidB12() {// digit with lowest count of unsolved
	int mincount =100, digmin = 10;
	for (uint32_t idig = 0; idig < 9; idig++) {
		int cc = FD[idig].Count();
		if (cc < 7) continue;// digit solved
		if (cc < mincount) {mincount = cc;	digmin = idig;	}
	}
	if (digmin > 8) return;// should never be unsolved status
	BF64 mysol = zh2b_g.fd_sols[0][digmin];
	if ((FD[digmin] & mysol) == mysol) { // if all true possible try it
		(this + 1)->GuessGo(digmin, mysol);
		if (zh2b_g.ua_ret) return;// return if  ua found
	}
	// try false 
	int r_unsolved = (rows_unsolved.bf.u64 >> zh2b_t_runsolvedshift[digmin]) & 077;
	BF64 tuaw[500], tsolw[500];
	int nuaw = zh2b1d_g.Go(mysol, FD[digmin], tsolw, tuaw, r_unsolved);
	//cout  << " nfalse perms=" << nuaw << endl;
	for (int i = 0; i < nuaw ; i++) {
		(this + 1)->GuessGo(digmin, tsolw[i]);
		if (zh2b_g.ua_ret) return;// return if  ua found
	}
}


void ZH2B::Debug(int all) {
	cout << "DEBUG  nbsol=" << zh2b_g.nsol << " unsolved=" << Unsolved_Count() << endl;
	//	cout << zh1b_g.out27 << " band1 "<<endl;
	char zi[82];  SetKnown(zi);
	cout << zi << " known rows 1_6 digits " << endl;
	if (!all) return;

	cout << "map per digit bands 2 3" << endl;
	for (int ib = 0; ib < 2; ib++) {
		for (int ir = 0; ir < 3; ir++) {
			for (int idig = 0; idig < 9; idig++) {
				unsigned vf = FD[idig].bf.u32[ib];
				unsigned wr = (vf >> (9 * ir)) & 0x1ff;
				for (int k = 0; k < 9; k++) {
					if (wr & (1 << k))		cout << idig + 1;
					else 		cout << ".";
					if (k == 2 || k == 5) 	cout << " ";
				}
				cout << "  ";
			}
			cout << endl; //end of row
		}
		cout << endl; // end of block
	}
	cout << endl; // end of map per digit

}
int ZH2B::GetAllDigits(int cell) {
	int ir = 0, xcell = C_To128[cell];;
	for (int i = 0; i < 9; i++) if (FD[i].On(xcell))ir |= (1 << i);
	return ir;
}
void ZH2B::ImageCandidats() {
	int dig_cells[81]; for (int i = 0; i < 54; i++) dig_cells[i] = GetAllDigits(i);
	int i, j, l, lcol[9], tcol = 0, ncand = 0;
	cout << "PM map " << endl << endl;
	for (i = 0; i < 9; i++) {  // attention ici i indice colonne
		lcol[i] = 2;    // 2  mini tous chiffres imposés
		for (j = 0; j < 6; j++) {
			l = _popcnt32(dig_cells[9 * j + i]);
			if (l > lcol[i])       lcol[i] = l;
		}
		tcol += lcol[i];
	}
	for (i = 0; i < 9; i++) {
		if ((i == 3) || (i == 6))cout << "|";
		cout << (char)('A' + i) << Blancs(lcol[i], 1);
	}
	cout << endl;
	for (i = 0; i < 6; i++) { // maintenant indice ligne
		if ((i == 3) || (i == 6)) {
			for (int ix = 0; ix < (tcol + 10); ix++)       cout << (char)'-';
			cout << endl;
		}
		for (j = 0; j < 9; j++) {
			if ((j == 3) || (j == 6))cout << "|";
			int cell = 9 * i + j, digs = dig_cells[cell], 
				ndigs = _popcnt32(digs);
			ncand += ndigs;
			for (int id = 0; id < 9; id++)if (digs & (1 << id))
				cout << id + 1;
			cout << Blancs(lcol[j] + 1 - ndigs, 1);
		} // end for j
		cout << endl;
	} // end for i
	cout << endl;

}



//============================= ZH2B5 uas guas creation

/*
process using as source
	fd_sols[1] is pm if source=0  fd_revised if  source=1
the solution is known and is fd_sols[0]
all digits not in bitfield fl are known (usualy fl 2-4 digits)
internal myfd as initial value  cleaned in zb2b_1d step by step

*/
uint64_t  ZH2B5_GLOBAL::FindUAsInit(int fl, int source) {
	BF64 * mypm = zh2b_g.fd_sols[1];
	if (source) mypm = zh2b_g.fd_revised;
	uint64_t solved_cells = 0;
	uint32_t nd = 0;
	for (int idig = 0, bit = 1; idig < 9; idig++, bit <<= 1) {
		if (fl & bit) {
			fdsw[0][nd] = zh2b_g.fd_sols[0][idig];
			fdsw[2][nd].bf.u64 = (~fdsw[0][nd].bf.u64) & BIT_SET_2X;
			myfd[nd++] = mypm[idig];
		}
		else solved_cells |= zh2b_g.fd_sols[0][idig].bf.u64;
	}
	ndigits = nd;
	cells_unsolved.bf.u64 = BIT_SET_2X ^ solved_cells;
	for (uint32_t i = 0; i < nd; i++) {
		myfd[i] &= cells_unsolved;
		if (myfd[i].Count() == 6)return 0;
	}
	return solved_cells;
	// first cleaning return unsolved cells in bits 	
}
void ZH2B5_GLOBAL::CollectUas5() {
	if (diag)cout << "entry CollectUas5()" << endl;
	nuaf5 = 0;
	memset(&zh2b5[0], 0, sizeof zh2b5[0]);
	memcpy(zh2b5[0].FD, myfd, sizeof myfd);
	memset(zh2b5[0].CompFD, 0, sizeof myfd);
	zh2b5[0].cells_unsolved = cells_unsolved;
	uint32_t t[5] = { 077,07777,0777777,077777777,07777777777 };
	zh2b5[0].rows_unsolved.bf.u32[0] = t[ndigits-1];
	if(diag)zh2b5[0].ImageCandidats();//<<<<<<<<<<<<<<<<<<<<<<<<<<
	//if (1) return;
	zh2b5[0].ComputeNext5();
}
void ZH2B5_GLOBAL::ValidSol5(uint64_t * sol) {//Ua to catch
	if (zh2b5_g.diag) {
		cout << "valid sol nuaf5="<<nuaf5<< endl;
	}
	uint64_t ua = 0, *s0 = &fdsw[0][0].bf.u64;
	for (uint32_t i = 0; i < ndigits; i++) {
		uint64_t w = sol[i] & ~s0[i];// digit ua
		if (!w) {
			if (zh2b5_g.diag) cout << "subset assumed rdigit" << i << endl;
			return; //a subset exists
		}
		ua |= w;
	}
	/* not true can be small gua not a band ua
	if (!modevalid) {// not gua mode must be 1n both bands
		register uint64_t R = BIT_SET_27;
		if (!(ua&R))return;; // ua in band2
		R <<= 32;
		if (!(ua&R))return; // ua in band1
	}
	*/
	uint64_t  cc = (int)_popcnt64(ua);
	if (zh2b5_g.diag)cout << "cc=" << cc << endl;
	if (cc > sizef5) return;
	if (nuaf5 >= 40 && modevalid) {// reasonnable limit for a given step
		int limit = tuaf5[39].bf.u64 >> 59;
		if (cc >= limit)return; // too many uas found here, skip it
		nuaf5 = 39;
	}
	ua |= cc << 59;
	//int ir = genuasb12.AddUA64(&tuaf5[0].bf.u64, nuaf5);
	//uint64_t cc = _popcnt64(genuasb12.ua);
	//genuasb12.ua.bf.u64  |= cc << 59;
	if (zh2b5_g.diag)cout << Char2Xout(ua) << " ua found nuaf5=" << nuaf5 << "cc=" << cc << endl;
	AddUA64(&tuaf5[0].bf.u64, nuaf5,ua);
}

//======================== ZH2B5  2-5 digits
int ZH2B5::Update5() {
	if(zh2b5_g.ndigits>5)return 0; // force false if bad use
	int Shrink = 1;
	register int S, A;
	register unsigned int cl, *wcl = FD[0].bf.u32;
	while (Shrink) {
		Shrink = 0;
		if (!rows_unsolved.bf.u32[0])return 1;// solved

		{register unsigned int  AR = rows_unsolved.bf.u32[0];// valid for digits 0,1,2,3,4
		if (!(AR & 077))goto digit1;

		//=digit 0
		if (FD[0].bf.u32[0] == CompFD[0].bf.u32[0])goto digit0b;
		UPDN(0, 0)	if ((AR & 7) != S) {
			AR &= 07777777770 | S;	UPWCL5(0, 2, 4, 6, 8)
		}

	digit0b:if (FD[0].bf.u32[1] == CompFD[0].bf.u32[1])goto digit1;
		UPDN(0, 1)	if (((AR >> 3) & 7) != S) {
			AR &= 07777777707 | (S << 3);	UPWCL5(1, 3, 5, 7, 9)
		}

	digit1:	if (!(AR & 07700))goto digit2;

		if (FD[1].bf.u32[0] == CompFD[1].bf.u32[0])goto digit1b;
		UPDN(1, 0)	if (((AR >> 6) & 7) != S) {
			AR &= 07777777077 | (S << 6); UPWCL5(0, 0, 4, 6, 8)
		}

	digit1b:if (FD[1].bf.u32[1] == CompFD[1].bf.u32[1])goto digit2;
		UPDN(1, 1)		if (((AR >> 9) & 7) != S) {
			AR &= 07777770777 | (S << 9); UPWCL5(1, 1, 5, 7, 9)
		}

	digit2:	
		if (zh2b5_g.ndigits < 3) goto end01234;
		if (!(AR & 0770000))goto digit3;

		if (FD[2].bf.u32[0] == CompFD[2].bf.u32[0])goto digit2b;
		UPDN(2, 0)	if (((AR >> 12) & 7) != S) {
			AR &= 07777707777 | (S << 12);	UPWCL5(0, 0, 2, 6, 8)
		}

	digit2b:if (FD[2].bf.u32[1] == CompFD[2].bf.u32[1])goto digit3;
		UPDN(2, 1)	if (((AR >> 15) & 7) != S) {
			AR &= 07777077777 | (S << 15);	UPWCL5(1, 1, 3, 7, 9)
		}

	digit3: 
		if (zh2b5_g.ndigits < 4) goto end01234;
		if (!(AR & 077000000))goto digit4;

		if (FD[3].bf.u32[0] == CompFD[3].bf.u32[0])goto digit3b;
		UPDN(3, 0)	  if (((AR >> 18) & 7) != S) {
			AR &= 07770777777 | (S << 18); UPWCL5(0, 0, 2, 4, 8)
		}

	digit3b:  if (FD[3].bf.u32[1] == CompFD[3].bf.u32[1])goto digit4;
		UPDN(3, 1)if (((AR >> 21) & 7) != S) {
			AR &= 07707777777 | (S << 21); UPWCL5(1, 1, 3, 5, 9)
		}

	digit4:if (zh2b5_g.ndigits < 5) goto end01234;
		if (!(AR & 07700000000))goto end01234;

		if (FD[4].bf.u32[0] == CompFD[4].bf.u32[0])goto digit4b;
		UPDN(4, 0)if (((AR >> 24) & 7) != S) {
			AR &= 07077777777 | (S << 24);  UPWCL5(0, 0, 2, 4, 6)
		}

	digit4b:if (FD[4].bf.u32[1] == CompFD[4].bf.u32[1])goto end01234;
		UPDN(4, 1)if (((AR >> 27) & 7) != S) {
			AR &= 0777777777 | (S << 27);  UPWCL5(1, 1, 3, 5, 7)
		}

	end01234: rows_unsolved.bf.u32[0] = AR;
		}// end of validity for AR 01234

	}// end while
	return 1;
}
int ZH2B5::FullUpdate5() {
	while (1) {
		if (!Update5()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (zh2b5_g.ndigits > 3) {
			if (ApplySingleOrEmptyCells5())	return 0; // locked 
			if (zh2b5_g.single_applied) continue;
		}
		break;
	}
	return 1;
}
int ZH2B5::ApplySingleOrEmptyCells5() {
#define NAKED5(X) Map=map[X];R2|=R1&Map;R1|=Map 
	zh2b5_g.single_applied = 0;
	uint64_t * map = &FD[0].bf.u64;
	uint64_t unsolved = cells_unsolved.bf.u64;
	register uint64_t R2 = map[0] & map[1],
		R1 = (map[0] | map[1]), Map;// digits 12
	if (zh2b5_g.ndigits > 2)NAKED5(2);
	if (zh2b5_g.ndigits > 3)NAKED5(3);
	if (zh2b5_g.ndigits > 4)NAKED5(4);
	if (unsolved & (~R1)) return 1; // locked
	R1 &= ~R2;
	R1 &= unsolved; // forget solved seen as singles
	if (!R1) return 0;
	zh2b5_g.single_applied = 1;

	while (R1) {// usually a very small number of cells to assign
		uint32_t  res;
		if (!bitscanforward64(res, R1)) break;
		uint64_t bit = (uint64_t)1 << res; // switch to the bit value
		R1 &= ~bit;  // clear the bit
		// call Seta(int digit, int xcell) so find digit
		for (uint32_t idig = 0; idig < zh2b5_g.ndigits; idig++) {
			if (map[idig] & bit) {// this is the digit
				if (FD[idig].Off(res))  return 1; // invalid, gane locked
				Seta5(idig, res);
				goto nextr1;// finished for that cell
			}
		}
		return 1; //conflict with a previous cell assugn
	nextr1: {}
	}
	return 0;// not locked
}
int ZH2B5::Seta5(int digit, int xcell) { // single in cell
	int cell = From_128_To_81[xcell],
		block = TblBoard_Block[cell];
	if (FD[digit].Off(xcell)) return 1; // not valid
	//Assign(digit, cell, xcell);
	BF64 *Fd = &FD[digit];
	*Fd &= AssignMask_Digit[cell].u64[0];
	int ddig = 6 * digit;
	rows_unsolved.Clear(ddig + C_row[cell]);//6*digit + row
	cells_unsolved.Clear(xcell);
	BF64 * RF = &FD[zh2b5_g.ndigits - 1];//last used digit
	for (; RF >= FD; RF--)RF->Clear(xcell);
	Fd->Set(xcell); // restore bit for digit assigned
	return 0;
}
void  ZH2B5::Guess5() {// solve in table digit with lowest count
	if (zh2b5_g.diag) {
		cout << oct << rows_unsolved.bf.u64 << dec << " unsolved guess5" << endl;
		ImageCandidats();
	}
	if (!rows_unsolved.bf.u64) {// this is a solution
		//cout << "valid sol" << endl;
		zh2b5_g.ValidSol5(&FD[0].bf.u64);
		return;
	}
	int mincount = 100, digmin,ncount=0;
	for (uint32_t idig = 0; idig < zh2b5_g.ndigits; idig++) {
		int cc = FD[idig].Count();
		if (cc < 7) continue;
		ncount++;
		if (cc < mincount) {
			mincount = cc;
			digmin = idig;
		}
	}
	if (zh2b5_g.diag) {
		cout << oct << rows_unsolved.bf.u64 << dec << " unsolved guess5 for dig"
			<<digmin+1<< endl;
		//ImageCandidats();
	}
	// put in table all digmin valid solutions 
	BF64 tuaw[100], tsolw[100];
	//cout << "call solve 1 digit for digitw=" << digmin+1 << endl;
	int nuaw=zh2b1d_g.Go(zh2b5_g.fdsw[0][digmin], FD[digmin], tsolw, tuaw,
		(rows_unsolved.bf.u32[0]) >> (6 * digmin) & 077);
	rows_unsolved.bf.u32[0] &= ~(077 << (6 * digmin));// clear unsolved
	if (zh2b5_g.diag)cout << "return nuaw=" << nuaw << " digmin="<<digmin+1<< endl;
	//try each digit full solution if not fully true (subset would exist)
	for (int i = 0; i < nuaw; i++) {
		ZH2B5 * nextz = this + 1;// open a new zh2b
		*nextz = *this;
		{// clear the solution for other digits
			register uint64_t R = tsolw[i].bf.u64, nR = ~R;
			if (zh2b5_g.diag)cout<<Char2Xout(R) << "apply for i=" << i << " digmin=" << digmin+1 << endl;
			uint64_t * RF = &nextz->FD[zh2b5_g.ndigits - 1].bf.u64;
			for (; RF >= &nextz->FD->bf.u64; RF--)*RF &= nR;
			nextz->FD[digmin].bf.u64 = R;// restore the solution for digit
		}	
		// if ncount=2 it is a solution
		if (ncount == 2) {
			if (zh2b5_g.diag) nextz->ImageCandidats();
			zh2b5_g.ValidSol5(&nextz->FD[0].bf.u64);
		}
		else nextz->ComputeNext5();// continue solving
	}
}
int ZH2B5::GetAllDigits(int cell) {
	int ir = 0, xcell = C_To128[cell];;
	for (uint32_t i = 0; i < zh2b5_g.ndigits; i++) if (FD[i].On(xcell))ir |= (1 << i);
	return ir;
}
void ZH2B5::ImageCandidats() {
	int dig_cells[81]; for (int i = 0; i < 54; i++) dig_cells[i] = GetAllDigits(i);
	int i, j, l, lcol[9], tcol = 0, ncand = 0;
	cout << "PM map " << endl << endl;
	for (i = 0; i < 9; i++) {  // attention ici i indice colonne
		lcol[i] = 2;    // 2  mini tous chiffres imposés
		for (j = 0; j < 6; j++) {
			l = _popcnt32(dig_cells[9 * j + i]);
			if (l > lcol[i])       lcol[i] = l;
		}
		tcol += lcol[i];
	}
	for (i = 0; i < 9; i++) {
		if ((i == 3) || (i == 6))cout << "|";
		cout << (char)('A' + i) << Blancs(lcol[i], 1);
	}
	cout << endl;
	for (i = 0; i < 6; i++) { // maintenant indice ligne
		if ((i == 3) || (i == 6)) {
			for (int ix = 0; ix < (tcol + 10); ix++)       cout << (char)'-';
			cout << endl;
		}
		for (j = 0; j < 9; j++) {
			if ((j == 3) || (j == 6))cout << "|";
			int cell = 9 * i + j, digs = dig_cells[cell],
				ndigs = _popcnt32(digs);
			ncand += ndigs;
			for (int id = 0; id < 9; id++)if (digs & (1 << id))
				cout << id + 1;
			cout << Blancs(lcol[j] + 1 - ndigs, 1);
		} // end for j
		cout << endl;
	} // end for i
	cout << endl;

}


//======================== ZH2B1B one digit 2 bands table of solutions

int ZH2B_1D_GLOBAL::Go(BF64 & sol, BF64 & fde, BF64 *tsol, BF64 *tua,int ru) {
	tsolw = tsol;
	tuaw = tua;
	mysol = sol;
	nsolw = 0;
	myandsol.bf.u64 = BIT_SET_2X;
	zh2b1d[0].FD = fde;
	zh2b1d[0].CompFD.bf.u64 = 0;
	zh2b1d[0].ComputeNext(ru);
	return nsolw;
}
int ZH2B_1D::GetSols( int ru) {
	CompFD.bf.u64 = 0;
	FD = zh2b_g.mystart &zh2b_g.myandsol;
	zh2b_g.nsolw = 0;
	ComputeNext(ru);
	return zh2b_g.nsolw;
}
int ZH2B_1D::GetAllSols(BF64 & fde, int ru, BF64 & fdsol) {
	zh2b_g.mystart = FD = fde;
	zh2b_g.mysol = fdsol;
	zh2b_g.myandsol.bf.u64 = BIT_SET_2X;
	zh2b_g.nsolw = 0;
	ComputeNext(ru);
	return zh2b_g.nsolw;
}
void ZH2B_1D::ComputeNext(int ru) {
	if (Update(ru)) {
		if (!ru) {// new sol put it in table
			if (FD != zh2b1d_g.mysol) {
				zh2b1d_g.tuaw[zh2b1d_g.nsolw] = FD- zh2b1d_g.mysol;// partial ua
				zh2b1d_g.tsolw[zh2b1d_g.nsolw++] = FD;// partial solution
			}
		}
		else Guess(ru);
	}
}
int ZH2B_1D::Update(int &ru) {// solve as much as possible 
	int Shrink = 1;
	register int S, A;
	while (Shrink) {
		Shrink = 0;
		if (!(ru & 7)) goto band2;
		A = FD.bf.u32[0];
		if (A == CompFD.bf.u32[0]) goto band2;
		Shrink = (TblShrinkMask[A & 0x1FF] |
			TblShrinkMask[(A >> 9) & 0x1FF] << 3 |
			TblShrinkMask[(A >> 18) & 0x1FF] << 6);
		if ((A &= TblComplexMask[Shrink]) == 0)  return 0;
		S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);
		FD.bf.u32[1] &= TblMaskSingle[S];
		S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];
		CompFD.bf.u32[0] = 
			FD.bf.u32[0] = A;
		ru = (ru & 070) | S;
	band2:;
		if (!(ru & 070)) goto exit;
		A = FD.bf.u32[1];
		if (A == CompFD.bf.u32[1]) goto exit;
		Shrink = (TblShrinkMask[A & 0x1FF] |
			TblShrinkMask[(A >> 9) & 0x1FF] << 3 |
			TblShrinkMask[(A >> 18) & 0x1FF] << 6);
		if ((A &= TblComplexMask[Shrink]) == 0)  return 0;
		S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);
		FD.bf.u32[0] &= TblMaskSingle[S];
		S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];
		CompFD.bf.u32[1] = 
			FD.bf.u32[1] = A;
		ru = (ru & 7) | (S<<3);
	exit:;
	}
	//cout << Char2Xout(FD.bf.u64) << " sortie update ru=" 
	//	<< oct << ru << dec << endl;

	return 1;
}
void ZH2B_1D::Guess(int ru) {//
	int dcell = 0,v,ruw=ru;
	{
		register uint32_t a = FD.bf.u32[0], b = FD.bf.u32[1];
		if (!(ru & 7)) goto guess_b2;
		if ((ru & 070) && (_popcnt32(b) < _popcnt32(a))) goto guess_b2;
		// guess in band1
		ruw = (-ru) &ru;
		v = FD.bf.u32[0] & TblRowUnsolved[ruw];// unknown in last unknown row
		goto gogo;
	guess_b2:	// guess in band2
		dcell = 27;
		ruw=ru >> 3;
		ruw = (-ruw) &ruw;
		v = FD.bf.u32[1] & TblRowUnsolved[ruw];// unknown in last unknown row
		ruw <<= 3; //relocate bit to the right place
	}
gogo:;
	uint32_t cell;
	ru ^= ruw;// clear the row for next steps
	while (bitscanforward(cell, v)) {
		v ^= 1 << cell;
		cell += dcell;
		ZH2B_1D * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Assign(cell);
		mynext->ComputeNext(ru);
	}
}
int ZH2B_1D::IsValid(uint64_t v) {
	CompFD.bf.u64 = 0;
	FD.bf.u64 = v;
	int ru = 077;
	return Update(ru);
}

//=================================== ZHONE

//===================== now ZHone code 

ZHONE_GLOBAL::ZHONE_GLOBAL() {
	zsol  = 0; // no solution unless required buy the user
}
void ZHONE_GLOBAL::GetBand(int * b, uint32_t * t) {
	band0 = b; tua = t;  nua = 0;
	// build sol per digit and pm per digit at start
	memset(fd_sols, 0, sizeof fd_sols);
	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 3; j++) {
			int cell = 9 * j + i, dig = b[cell];
			fd_sols[1][dig] |= Zhoucol << i; // add candidates in the column
			fd_sols[0][dig] |= 1 << cell;
		}
	}
}
void ZHONE_GLOBAL::GetBand(uint32_t fd_solsb[2][9], int * b, uint32_t * t) {
	band0 = b; tua = t;  nua = 0;
	memcpy(fd_sols, fd_solsb, sizeof fd_sols);
}


void ZHONE_GLOBAL::AddUA(uint32_t ua) {// add if and clear supersets
	// ua 27 bit + 5 bit length
	if (nua >= 100) return; // size given in STD_B416 go_17_bands
	//if(diag)cout << "add in tua " << tua << " nua="<< nua<< endl;
	register uint32_t ua27 = ua & BIT_SET_27;
	for (uint32_t iua = 0; iua < nua; iua++) {
		register uint32_t R = tua[iua];
		if (R < ua) {// is it subset
			R&= BIT_SET_27;
			if ((R&ua27) == R) {
				//if (diag)cout << Char27out(R) << " subset of (iua=)"<<iua << endl
				//	<< Char27out(ua27) << endl;
				return;// we have a subset
			}
		}
		else if (R == ua)return;
		else {
			for (uint32_t jua = nua; jua > iua; jua--)tua[jua] = tua[jua - 1];// to insert the new
			tua[iua] = ua;// inserted
			nua++;
			//if (diag)cout << Char27out(ua27) << " inserted nua=" << nua << endl;
			// is it a subset of a previous entry
			for (iua++; iua < nua; iua++) {
				if ((tua[iua] &ua27) == ua27) {// we have a subset
					//if (diag)cout << Char27out(tua[iua]) << " superset cleared nua=" << nua << endl;
					for (uint32_t k = iua + 1; k < nua; k++)tua[k - 1] = tua[k];
					nua--;
					iua--; //continue same position
				}
			}
			return;
		}
	}
	tua[nua++] = ua;// added
	//if (diag) {
		//cout << Char27out(ua27) << " added nua=" << nua << endl;
		//cout << Char27out(tua[nua-1]) << " added stored"  << endl;
	//}
}
void ZHONE_GLOBAL::FindMissingUAs() {
	struct SPB3 {// spots to find band 3 missing uas
		int  possible_cells, all_previous_cells, active_cells, iuab3;
	}spb3[12], *s3, *sn3;
	s3 = spb3;
	memset(s3, 0, sizeof spb3[0]);
	s3->active_cells = BIT_SET_27;// all cells active
	s3->possible_cells = tua[0]& BIT_SET_27;
	int tcells[10], ispot;
	cout << "search nua=" << nua << endl;
	//____________________  here start the search
next:
	uint32_t cell;
	ispot = (int)(s3 - spb3); // be sure to have always the right one
	if (ispot > 10) return;
	if (!bitscanforward(cell, s3->possible_cells))goto back;
	{// apply cell in bitfields
		register int bit = (uint64_t)1 << cell;
		tcells[ispot] = cell;
		s3->possible_cells ^= bit;// clear bit
		s3->active_cells ^= bit;
		sn3 = s3 + 1;
		*sn3 = *s3;
		sn3->all_previous_cells |= bit;
		register int filter = sn3->all_previous_cells,
			ac = s3->active_cells;
		// nextspot:take the next available ua to loop
		for (uint32_t i = s3->iuab3 + 1; i < nua; i++) {
			if (tua[i] & filter)continue;
			//cout << Char27out(tua[i]) << "next ua i=" <<i << endl;
			sn3->iuab3 = i;
			sn3->possible_cells = tua[i] & ac;
			s3 = sn3; // switch to next spot
			goto next;
		}
		// no more ua is it valid or take more uas
		zhone[0].AddMissingUAs(tcells, ispot + 1);
		goto next;
	}
	// going back, for a non empty index, count it back
back:
	//if (1) return;
	if (--s3 >= spb3)goto next;

}
void ZHONE_GLOBAL::ValidPuzzle(uint32_t * sol) {//depending on type
	//cout << "entry valid puzzle" << endl;
	switch (type) {
	case 1: {// mode add uas
		uint32_t or_sol = 0;
		//for (int i = 0; i < 9; i++) {// 9 digits
		for (int i = 0; i <9; i++) {// 9 digits
			or_sol |= sol[i] & ~(fd_sols[0][i]);// invalid cells
		}
		if (or_sol) {
			//cout << Char27out(or_sol) << " ua found "<< _popcnt32(or_sol) << endl;
			AddUA(or_sol | _popcnt32(or_sol) << 27);
		}
	}
		break;
	}
	return;
}

void ZHONE_GLOBAL::PrintTua() {
	cout << "band status" << endl;
	for (uint32_t i = 0; i < nua; i++) {
		uint32_t ua = tua[i], ua27 = ua & BIT_SET_27, nd = ua >> 27;
		cout << i << "\t" << Char27out(ua27) << "\t" << nd << endl;
	}
}


//============================================ ZHONE
void ZHONE::InitOne_std_band() {//init after zh1b_g getband
	cells_unsolved = BIT_SET_27;
	memcpy(FD, zh1b_g.fd_sols[1], sizeof FD);
	for (int i = 0; i < 9; i++)	FD[i] |= 7 << 27;//set unknown rows
	zh1b_g.ndigits = 9;
}
void ZHONE::CheckSolPerDigit() {//init after zh1b_g getband
	cells_unsolved = BIT_SET_27;
	memcpy(FD, zh1b_g.fd_sols[0], sizeof FD);
	zh1b_g.ndigits = 9;
	ImageCandidats();
}


int ZHONE::InitSudokux(GINT * t, int n) {// cells must be 0-26
	int Digit_cell_Assigned[9];
	memset(Digit_cell_Assigned, 0, sizeof Digit_cell_Assigned);
	*this=zhone_i;
	for (int icell = 0; icell < n; icell++) {
		int cell = t[icell].u8[0], digit = t[icell].u8[1], bit = 1 << cell;
		if ((FD[digit] & bit) == 0)  return 1;// check not valid entry
		Assign(digit, cell);
		Digit_cell_Assigned[digit] |= bit;
	}
	for (int i = 0; i < 9; i++)  FD[i] &=
		cells_unsolved |
		Digit_cell_Assigned[i] |
		07000000000;// don't touch rows unsolved

//	Debug(1);
//	ImageCandidats();
	return 0;
}
int ZHONE::Isvalid() { // usually after init 2 steps
	zh1b_g.InitIsvalid();
	ComputeNext();
	return zh1b_g.nsol;
}

#define NAKED(X) 	Map=map[X];R3|=R2&Map;R2|=R1&Map;R1|=Map;
int ZHONE::ApplySingleOrEmptyCells() {
	zh1b_g.single_applied = 0;
	uint32_t * map = FD, unsolved = cells_unsolved;
	register int R2 = map[0] & map[1],
		R1 = (map[0] | map[1]), Map, R3;// digits 12
	Map = map[2]; R3 = R2 & Map; R2 |= R1 & Map; R1 |= Map;
	NAKED(3) NAKED(4) NAKED(5) NAKED(6)	NAKED(7) NAKED(8) // digits 3-9
		if (unsolved & (~R1)) return 1; // locked
	R1 &= ~R2;
	R1 &= unsolved; // forget solved seen as singles
	if (R1) zh1b_g.single_applied = 1;
	else {
		zh1b_g.pairs = R2 & (~R3);
		zh1b_g.triplets = R3;
		return 0;
	}
	while (R1) {// usually a very small number of cells to assign
		uint32_t res;
		if (!bitscanforward(res, R1)) break;
		int cell = res, bit = 1 << cell; // switch to the bit value
		//		cout << "naked cell" << From_128_To_81[res]<< endl;
		R1 &= ~bit;  // clear the bit
		// call Seta(int digit, int xcell) so find digit 
		for (int idig = 0; idig < 9; idig++) {
			if (map[idig] & bit) {// this is the digit
				//				if (FD[idig].Off(res))  return 1; // invalid, gane locked
				//				Seta(idig, res);
				int cell = res;
				Assign(idig, cell);
				goto nextr1;// finished for that cell
			}
		}
		return 1; //conflict with a previous cell assign
	nextr1: {}
	}
	return 0;// not locked 
}
void ZHONE::Seta(int digit, int cell) { // single in cell
	Assign(digit, cell);
	register  int  bit= 1 << cell,cl = ~bit;
	FD[0]&=cl; FD[1] &= cl; FD[2] &= cl; FD[3] &= cl; FD[4] &= cl;
	FD[5] &= cl; FD[6] &= cl; FD[7] &= cl; FD[8] &= cl;
	FD[digit] |= bit;// restore digit
}
int ZHONE::Update() {
	int Shrink = 1;
	register int S, A, cl;
	register uint32_t *wcl = FD;
	while (Shrink) {
		Shrink = 0;
		if ((A = FD[0]) - CompFD[0]) { UPDN1(0, 1, 2, 3, 4, 5, 6, 7, 8); }
		if ((A = FD[1]) - CompFD[1]) { UPDN1(1, 0, 2, 3, 4, 5, 6, 7, 8); }
		if ((A = FD[2]) - CompFD[2]) { UPDN1(2, 0, 1, 3, 4, 5, 6, 7, 8); }
		if ((A = FD[3]) - CompFD[3]) { UPDN1(3, 0, 1, 2, 4, 5, 6, 7, 8); }
		if ((A = FD[4]) - CompFD[4]) { UPDN1(4, 0, 1, 2, 3, 5, 6, 7, 8); }
		if ((A = FD[5]) - CompFD[5]) { UPDN1(5, 0, 1, 2, 3, 4, 6, 7, 8); }
		if ((A = FD[6]) - CompFD[6]) { UPDN1(6, 0, 1, 2, 3, 4, 5, 7, 8); }
		if ((A = FD[7]) - CompFD[7]) { UPDN1(7, 0, 1, 2, 3, 4, 5, 6, 8); }
		if ((A = FD[8]) - CompFD[8]) { UPDN1(8, 0, 1, 2, 3, 4, 5, 6, 7); }
		//	  Debug(1);
	}// end while
#ifdef DIAG
	cout << "end update cycle" << endl;
	Debug(1);
	ImageCandidats();
#endif
	return 1;
}
int ZHONE::Update4() {
	int Shrink = 1;
	register int S, A, cl;
	register uint32_t *wcl = FD;
	while (Shrink) {
		Shrink = 0;
		if ((A = FD[4]) - CompFD[4]) { UPDN1(4, 0, 1, 2, 3, 5, 6, 7, 8); }
		if ((A = FD[3]) - CompFD[3]) { UPDN1(3, 0, 1, 2, 4, 5, 6, 7, 8); }
		if ((A = FD[2]) - CompFD[2]) { UPDN1(2, 0, 1, 3, 4, 5, 6, 7, 8); }
		if ((A = FD[1]) - CompFD[1]) { UPDN1(1, 0, 2, 3, 4, 5, 6, 7, 8); }
		if ((A = FD[0]) - CompFD[0]) { UPDN1(0, 1, 2, 3, 4, 5, 6, 7, 8); }
	}// end while
	return 1;
}
int ZHONE::Update6() {
	int Shrink = 1;
	register int S, A, cl;
	register uint32_t *wcl = FD;
	while (Shrink) {
		Shrink = 0;
		if ((A = FD[5]) - CompFD[5]) { UPDN1(5, 0, 1, 2, 3, 4, 6, 7, 8); }
		if ((A = FD[4]) - CompFD[4]) { UPDN1(4, 0, 1, 2, 3, 5, 6, 7, 8); }
		if ((A = FD[3]) - CompFD[3]) { UPDN1(3, 0, 1, 2, 4, 5, 6, 7, 8); }
		if ((A = FD[2]) - CompFD[2]) { UPDN1(2, 0, 1, 3, 4, 5, 6, 7, 8); }
		if ((A = FD[1]) - CompFD[1]) { UPDN1(1, 0, 2, 3, 4, 5, 6, 7, 8); }
		if ((A = FD[0]) - CompFD[0]) { UPDN1(0, 1, 2, 3, 4, 5, 6, 7, 8); }
	}// end while
	return 1;
}
int ZHONE::Update7() {
	int Shrink = 1;
	register int S, A, cl;
	register uint32_t *wcl = FD;
	while (Shrink) {
		Shrink = 0;
		if ((A = FD[6]) - CompFD[6]) { UPDN1(6, 0, 1, 2, 3, 4, 5, 7, 8); }
		if ((A = FD[5]) - CompFD[5]) { UPDN1(5, 0, 1, 2, 3, 4, 6, 7, 8); }
		if ((A = FD[4]) - CompFD[4]) { UPDN1(4, 0, 1, 2, 3, 5, 6, 7, 8); }
		if ((A = FD[3]) - CompFD[3]) { UPDN1(3, 0, 1, 2, 4, 5, 6, 7, 8); }
		if ((A = FD[2]) - CompFD[2]) { UPDN1(2, 0, 1, 3, 4, 5, 6, 7, 8); }
		if ((A = FD[1]) - CompFD[1]) { UPDN1(1, 0, 2, 3, 4, 5, 6, 7, 8); }
		if ((A = FD[0]) - CompFD[0]) { UPDN1(0, 1, 2, 3, 4, 5, 6, 7, 8); }
	}// end while
	return 1;
}
void ZHONE::Guess() {
	if (!cells_unsolved) {
		if (zh1b_g.type) {
			//ImageCandidats();
			zh1b_g.ValidPuzzle(FD); return;	}
		else if (zh1b_g.zsol && (!zh1b_g.nsol)) SetKnown(zh1b_g.zsol);// store the first solution
		zh1b_g.nsol++;
		if (zh1b_g.nsol > zh1b_g.lim) zh1b_g.go_back = 1;
		return;
	}
	uint32_t res;
	int ndig = 2;
	register int R3 = zh1b_g.pairs;
	if (bitscanforward(res, R3))goto gogo; // first pair is ok to go
	ndig = 3;
	R3 = zh1b_g.triplets;
	bitscanforward(res, R3); // first triplet is ok to go
gogo:
	int cell = res, bit = 1 << cell;
	//cout << "brute force guess  cell " << cellsFixedData[cell].pt << endl;
	for (int idig = 0; idig < 9; idig++) {
		if (FD[idig] & bit) {// one valid digit		
			if (--ndig) {
				ZHONE * mynext = this + 1; // start next guess
				*mynext=*this;
				mynext->Seta(idig, res);
				mynext->ComputeNext();
				if (zh1b_g.go_back) return;
			}
			else {// this is the last do it in the same place
				Seta(idig, res);
				ComputeNext();
				return;
			}
		}
	}
}
int ZHONE::FullUpdate() {
	if (zh1b_g.go_back) return 0;
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (ApplySingleOrEmptyCells())			return 0; // locked empty cell or conflict singles in cells
		if (zh1b_g.single_applied) {
			//			cout << "after singles" << endl; Debug();
			continue;
		}
		break;
	}
	return 1;
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
char * ZHONE::SetKnown(char * zs) {
	strcpy(zs, &empty_puzzle[27]);
	zs[27] = 0;
	for (uint32_t idig = 0; idig < zh1b_g.ndigits; idig++) {// one digit
		int arows = FD[idig] >> 27;// 3 rows
		if (arows == 7) continue;
		int band = FD[idig];
		for (int j = 0; j < 3; j++) if (!(arows & (1 << j))) {
			int row = (band >> TblMult9[j]) & 0x1ff;
			uint32_t  irow;
			bitscanforward(irow, row);
			int	cell = TblMult9[j] + irow;
			zs[cell] = idig + '1';
		}
	}
	return zs;
}
void ZHONE::Debug(int all) {
	cout << "DEBUG  nbsol=" << zh1b_g.nsol << " unsolved=" << Unsolved_Count() 
		<< " index to zhone="<<this-zhone<< endl;
	//cout << zh1b_g.out27 << " band1 " << endl;
	char zi[82];  SetKnown(zi);
	cout << zi << " known rows 1_3 digits " << endl;
	if (!all) return;

	cout << "map per digit one band" << endl;
	for (int ir = 0; ir < 3; ir++) {
		for (uint32_t idig = 0; idig < zh1b_g.ndigits; idig++) {
			int vf = FD[idig], wr = (vf >> (9 * ir)) & 0x1ff;
			for (int k = 0; k < 9; k++) {
				if (wr & (1 << k))		cout << idig + 1;
				else 		cout << ".";
				if (k == 2 || k == 5) 	cout << " ";
			}
			cout << "  ";
		}
		cout << endl; //end of row
	}
	cout << endl; // end of map per digit

}
void ZHONE::DebugDigit(int digit) {
	cout << "DEBUG  digit=" << digit + 1 << endl;
	for (int ir = 0; ir < 3; ir++) {
		int vf = FD[digit], wr = (vf >> (9 * ir)) & 0x1ff;
		for (int k = 0; k < 9; k++) {
			if (wr & (1 << k))		cout << digit + 1;
			else 		cout << ".";
			if (k == 2 || k == 5) 	cout << " ";
		}
		cout << endl; //end of row
	}
	cout << endl; // end of block

}
int ZHONE::GetAllDigits(int cell) {
	int ir = 0, bit = 1 << cell;
	for (int i = 0; i < (int)zh1b_g.ndigits; i++) if (FD[i] & bit)ir |= (1 << i);
	return ir;
}
void ZHONE::ImageCandidats() {
	int dig_cells[27]; for (int i = 0; i < 27; i++) dig_cells[i] = GetAllDigits(i);
	int i, j, l, lcol[9], tcol = 0, ncand = 0;
	cout << "PM map " << endl << endl;
	for (i = 0; i < 9; i++) {  // column
		lcol[i] = 2;    // 2  mini 
		for (j = 0; j < 3; j++) {
			l = _popcnt32(dig_cells[9 * j + i]);
			if (l > lcol[i])       lcol[i] = l;
		}
		tcol += lcol[i];
	}
	for (i = 0; i < 9; i++) {
		if ((i == 3) || (i == 6))cout << "|";
		cout << (char)('A' + i) << Blancs(lcol[i], 1);
	}
	cout << endl;
	for (i = 0; i < 3; i++) { // now row 
		for (j = 0; j < 9; j++) {
			if ((j == 3) || (j == 6))cout << "|";
			int cell = 9 * i + j, digs = dig_cells[cell], 
				ndigs = _popcnt32(digs);
			ncand += ndigs;
			for (int id = 0; id < (int)zh1b_g.ndigits; id++)
				if (digs & (1 << id))
					if(zh1b_g.ndigits==9)	cout << id + 1;
					else cout << zh1b_g.digmap [id] + 1;
			cout << Blancs(lcol[j] + 1 - ndigs, 1);
		} // end for j
		cout << endl;
	} // end for i
	cout << endl;

}

//============== collector mode
void ZHONE::Checkstart() {
	zh1b_g.ndigits = 9;
	memcpy(FD, zh1b_g.fd_sols[1], sizeof FD);
	cells_unsolved = BIT_SET_27 ;
	for (uint32_t idig = 0; idig < zh1b_g.ndigits; idig++)
		FD[idig] |= 7 << 27;
	ImageCandidats();
}
int ZHONE::Start_nFloors(int floors) {
	uint32_t solved_cells = 0,nd=0;
	for (int idig = 0,bit=1; idig < 9; idig++,bit<<=1) {
		if (floors & bit) {
			zh1b_g.digmap[nd] = idig;
			zh1b_g.fdsw[0][nd]= zh1b_g.fd_sols[0][idig];
			zh1b_g.fdsw[2][nd] = (~zh1b_g.fdsw[0][nd]) & BIT_SET_27;
			FD[nd++] = zh1b_g.fd_sols[1][idig];
		}
		else solved_cells|= zh1b_g.fd_sols[0][idig];
	}
	zh1b_g.ndigits = nd;
	cells_unsolved=BIT_SET_27^ solved_cells;
	for (uint32_t i = 0; i < nd; i++) {
		FD[i] &= cells_unsolved;
		if (_popcnt32(FD[i]) == 3 && (!zh1b_g.modegua))return 1;
	}
	if (zh1b_g.diag) {
		cout << "end start nfloors" << endl;
		ImageCandidats();
	}
	return 0;
}
void ZHONE::Start_Uas_Mini(int floors, int floors_mini_row) {
	//cout << "start uas mini lfoors0" << oct << floors << " digits0" << floors_mini_row << dec << endl;
	zh1b_g.floors_mini_row = floors_mini_row; //used to check gangster changes 
	uint32_t solved_cells = 0, nd = 0;
	// first assigned digits
	for (int idig = 0, bit = 1; idig < 9; idig++, bit <<= 1) {
		if (!(floors & bit))  solved_cells |= zh1b_g.fd_sols[0][idig];
	}
	cells_unsolved = BIT_SET_27 ^ solved_cells;
	// then digits not exchanged
	floors &= ~floors_mini_row;
	for (int idig = 0, bit = 1; idig < 9; idig++, bit <<= 1) {
		if (floors & bit) {
			zh1b_g.digmap[nd] = idig; // used in gangster changes
			zh1b_g.digmap2[idig] = nd; // used in gangster changes
			zh1b_g.fdsw[0][nd] = zh1b_g.fd_sols[0][idig];
			zh1b_g.fdsw[2][nd] = (~zh1b_g.fdsw[0][nd]) & BIT_SET_27;
			FD[nd++] = zh1b_g.fd_sols[1][idig]& cells_unsolved;
		}
	}
	floors =floors_mini_row;// then digits  exchanged
	for (int idig = 0, bit = 1; idig < 9; idig++, bit <<= 1) {
		if (floors & bit) {
			zh1b_g.digmap[nd] = idig; // used in gangster changes
			zh1b_g.digmap2[idig] = nd; // used in gangster changes
			zh1b_g.fdsw[0][nd] = zh1b_g.fd_sols[0][idig];
			zh1b_g.fdsw[2][nd] = (~zh1b_g.fdsw[0][nd]) & BIT_SET_27;
			FD[nd++] = zh1b_g.fd_sols[1][idig] & cells_unsolved;
		}
	}
	zh1b_g.ndigits = nd;
	//cout << "end start uas nd=" << nd << endl;
	//ImageCandidats();
}
void ZHONE::ApplyGangsterChanges(int * g0, int * g1) {
	for (int i = 0,col=Zhoucol; i < 9; i++,col<<=1) {
		if (g0[i] == g1[i])continue;
		int changes = g0[i] ^ g1[i]; // one added one cleared
		// safety temp control, must be 2 digits exchanged
		if (changes & (~zh1b_g.floors_mini_row))return;
		//if (zh1b_g.diag > 1)
		for (int d = 0,bit=1; d < 9; d++,bit<<=1) {// check digits
			if (!(changes & bit)) continue;
			int digit = zh1b_g.digmap2[d];// digit rank in the brute force
			if (g0[i] & bit)FD[digit] &= ~col; 
			else FD[digit] |= col & cells_unsolved;
		}
	}
	if (zh1b_g.diag) {
		cout << "end apply gangster changes" << endl;
		ImageCandidats();
	}

}
void ZHONE::InitGuess() {// after start floors before guess
	for (uint32_t idig = 0; idig < zh1b_g.ndigits; idig++)
		FD[idig] |= 7 << 27;
	memset(zh1b_g.previous_ua_status, 0,
		sizeof zh1b_g.previous_ua_status); //no previous status
	for (int i = 0; i < 6; i++)
		zh1b_g.upstream_unsolved_cells[i] = BIT_SET_27;
	if (zh1b_g.diag) {
		cout << "end Init guess" << endl;
		ImageCandidats();
		cout << "end Init guess bis" << endl;
		ImageCandidats();
	}
}
int ZHONE::UpdateDigit(int digit) {
	register int S, A = FD[digit], cl;
	if (A  != CompFD[digit]) {
		int Shrink = (TblShrinkMask[A & 0x1FF] |
			TblShrinkMask[(A >> 9) & 0x1FF] << 3 |
			TblShrinkMask[(A >> 18) & 0x1FF] << 6);
		if ((A &= TblComplexMask[Shrink]) == 0)  return 0;
		S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);
		S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];
		if ((A >> 27) != S) {
			cl = ~(A & TblRowMask[S]);
			A = (A&BIT_SET_27) | (S << 27);
			cells_unsolved &= cl;
		}
		FD[digit] = CompFD[digit] = A;
	}// end if
	return 1;
}
void ZHONE::Set2(int cell) { // single in cell
	Assign(0, cell);
	register  int  cl = ~(1 << cell);
	FD[1] &= cl;
}
void ZHONE::Guess2() {	// note if all 1 are settled, the puzzle is solved
	if (zh1b_g.diag) {
		cout << Char27out(zh1b_g.previous_ua_status[0]) << "guess 2 with previous" << endl;
		cout << Char27out(zh1b_g.upstream_unsolved_cells[0]) << "upstream" << endl;
		ImageCandidats();	}
	int v = FD[0], vr = v >> 27;
	if (!vr) {
		register uint32_t w0 = FD[0] & zh1b_g.fdsw[2][0],
			w1 = FD[1] & zh1b_g.fdsw[2][1];
		if (!zh1b_g.modegua) 
			if ((!w0) || (!w1)) return;
		w0 |= (w1 | zh1b_g.previous_ua_status[0]);
		w0 |= _popcnt32(w0) << 27;
		zh1b_g.AddUA(w0);
		if (zh1b_g.diag)	cout << Char27out(w0) << " ua added direct guess2 + count " << (w0 >> 27) << "  nua=" << zh1b_g.nua << endl;
		return;
	}
	vr = (-vr)&vr; // catch last bit
	v &= TblRowUnsolved[vr];// unknown rows
	//cout << Char27out(v) << " v guess2" << endl;
	uint32_t cell;
	while (bitscanforward(cell, v)) {// first cell
		v ^= 1 << cell; // clear bit
		ZHONE * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Set2(cell);
		if (mynext->UpdateDigit(0)) {// update digit 1/0
			if (mynext->FD[0] >> 27) {
				mynext->Guess2();
			}
			else {// all digit 0 solved this is a solution
				//cout << Char27out(mynext->FD[0]) << " solution digit 0 " << endl;
				int unsolved = zh1b_g.upstream_unsolved_cells[0] & ~mynext->FD[0];
				register uint32_t w0 = mynext->FD[0] & zh1b_g.fdsw[2][0],
					w1 =( mynext->FD[1] & unsolved) & zh1b_g.fdsw[2][1];
				if ((!w0) || (!w1)) continue;
				w0 |= (w1  | 	zh1b_g.previous_ua_status[0]);
				w0 |= _popcnt32(w0) << 27;
				zh1b_g.AddUA(w0);
				if (zh1b_g.diag)	cout << Char27out(w0) << " ua added guess2 + count " << (w0 >> 27) << "  nua=" << zh1b_g.nua << endl;
			}
		}
	}
}
void ZHONE::Guess3() {	// solve digit 3/2
	if (zh1b_g.diag) {
		cout<<Char27out(zh1b_g.previous_ua_status[1]) << "guess 3 with previous" << endl;
		cout << Char27out(zh1b_g.upstream_unsolved_cells[1]) << "upstream unsolved" << endl;
		ImageCandidats();
	}
	int v = FD[2], vr = v >> 27;
	if (!vr) {
		uint32_t w3 = FD[2] & zh1b_g.fdsw[2][2];
		if ((!zh1b_g.modegua) && (!w3)) return;
		zh1b_g.upstream_unsolved_cells[0] =
			zh1b_g.upstream_unsolved_cells[1]  &~FD[2];
		if (zh1b_g.diag)			cout << Char27out(w3) << " guess3 ua" << endl;
		zh1b_g.previous_ua_status[0] = w3 | zh1b_g.previous_ua_status[1];
		//				mynext->FD[4] &= (mynext->cells_unsolved | 7 << 27);
		uint32_t cl = (zh1b_g.upstream_unsolved_cells[0] | 7 << 27);
		FD[1] &= cl; FD[0] &= cl;
		Guess2();
		return;
	}
	vr = (-vr)&vr; // catch last bit
	v &= TblRowUnsolved[vr];// unknown in last unknown row
	uint32_t cell;
	while (bitscanforward(cell, v)) {// first cell
		v ^= 1 << cell; // clear bit
		ZHONE * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Assign(2, cell);
		if (mynext->UpdateDigit(2)) {// update digit 3/2
			if (mynext->FD[2] >> 27) {
				mynext->Guess3();
			}
			else {// all digit 3 solved go to guess2 if not a zero deviation
				if (zh1b_g.diag > 1)
					cout << Char27out(mynext->FD[2]) << " guess2 " << endl;
				uint32_t w3 = mynext->FD[2] & zh1b_g.fdsw[2][2];
				if (!w3) continue;
				zh1b_g.upstream_unsolved_cells[0] =
					zh1b_g.upstream_unsolved_cells[1] & ~mynext->FD[2];
				zh1b_g.previous_ua_status[0] = w3 | zh1b_g.previous_ua_status[1];
				uint32_t cl = (zh1b_g.upstream_unsolved_cells[0] | 7 << 27);
				mynext->FD[1] &= cl;		mynext->FD[0] &= cl;
				mynext->Guess2();
			}
		}
	}
}
void ZHONE::Guess4() {	//  solve digit 4/3
	if (zh1b_g.diag) {
		cout << Char27out(zh1b_g.previous_ua_status[2]) << "guess 4" << endl;
		cout << Char27out(zh1b_g.upstream_unsolved_cells[2]) << "upstream unsolved" << endl;
		ImageCandidats();
	}
	int v = FD[3], vr = v >> 27;
	if (!vr) {
		uint32_t w4 = FD[3] & zh1b_g.fdsw[2][3];
		if ((!zh1b_g.modegua) && (!w4)) return;
		zh1b_g.upstream_unsolved_cells[1] =
			zh1b_g.upstream_unsolved_cells[2]  &~FD[3];
		if (zh1b_g.diag)			cout << Char27out(w4) << " direct guess4 ua" << endl;
		zh1b_g.previous_ua_status[1] = w4 | zh1b_g.previous_ua_status[2];
		//				mynext->FD[4] &= (mynext->cells_unsolved | 7 << 27);
		uint32_t cl = (zh1b_g.upstream_unsolved_cells[1] | 7 << 27);
		FD[2] &= cl; FD[1] &= cl; FD[0] &= cl;
		Guess3();
		return;
	}
	vr = (-vr)&vr; // catch last bit
	v &= TblRowUnsolved[vr];// unknown in last unknown row
	uint32_t cell;
	while (bitscanforward(cell, v)) {// first cell
		v ^= 1 << cell; // clear bit
		ZHONE * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Assign(3, cell);
		if (mynext->UpdateDigit(3)) {// update digit 4/3
			if (mynext->FD[3] >> 27) {
				mynext->Guess4();
			}
			else {// all digit 4 solved go to guess3 if not a zero deviation
				uint32_t w4 = mynext->FD[3] & zh1b_g.fdsw[2][3];
				if (!w4) continue;
				if (zh1b_g.diag > 1)
					cout << Char27out(mynext->FD[3]) << " new guess4" << endl;
				zh1b_g.upstream_unsolved_cells[1] =
					zh1b_g.upstream_unsolved_cells[2] & ~mynext->FD[3];
				zh1b_g.previous_ua_status[1] = w4 | zh1b_g.previous_ua_status[2];
				uint32_t cl = (zh1b_g.upstream_unsolved_cells[1] | 7 << 27);
				mynext->FD[2] &= cl;
				mynext->Guess3(); mynext->FD[1] &= cl;	mynext->FD[0] &= cl;
			}
		}
	}
}
void ZHONE::Guess5() {	//  solve digit 5/4
	if (zh1b_g.diag) {		
		cout << Char27out(zh1b_g.previous_ua_status[3]) << "guess 5 previous ua" << endl;
		cout << Char27out(zh1b_g.upstream_unsolved_cells[3]) << "upstream unsolved" << endl;
		ImageCandidats();	}
	if (!Update4()) return;;// digit 5 done in guess6
	int v = FD[4], vr = v >> 27;
	if (!vr) {
		uint32_t w5 = FD[4] & zh1b_g.fdsw[2][4];
		if ((!zh1b_g.modegua) && (!w5)) return;
		zh1b_g.upstream_unsolved_cells[2] =
			zh1b_g.upstream_unsolved_cells[3] & ~FD[4];
		if (zh1b_g.diag)			cout << Char27out(w5) << " direct guess5 ua" << endl;
		zh1b_g.previous_ua_status[2] = w5 | zh1b_g.previous_ua_status[3];
		uint32_t cl = (zh1b_g.upstream_unsolved_cells[2] | 7 << 27);
		FD[3] &= cl; FD[2] &= cl; FD[1] &= cl; FD[0] &= cl;
		Guess4();
		return;
	}
	vr = (-vr)&vr; // catch last bit
	v &= TblRowUnsolved[vr];// unknown in last unknown row
	uint32_t cell;
	while (bitscanforward(cell, v)) {// first cell
		v ^= 1 << cell; // clear bit
		ZHONE * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Assign(4, cell);
		if (mynext->UpdateDigit(4)) {// update digit 4/3
			if (mynext->FD[4] >> 27) {
				mynext->Guess5();
			}
			else {// all digit 5 solved go to guess4 if not a zero deviation
				uint32_t w5 = mynext->FD[4] & zh1b_g.fdsw[2][4];
				if (!w5) continue;
				if (zh1b_g.diag) {
					cout << Char27out(mynext->FD[4]) << " new guess5 " << endl;
					cout << Char27out(w5) << " new guess5 ua" << endl;
				}
				zh1b_g.upstream_unsolved_cells[2] =
					zh1b_g.upstream_unsolved_cells[3]  &~mynext->FD[4];
				zh1b_g.previous_ua_status[2] = w5 | zh1b_g.previous_ua_status[3];
				uint32_t cl = (zh1b_g.upstream_unsolved_cells[2] | 7 << 27);
				mynext->FD[3] &= cl;
				mynext->FD[2] &= cl;	mynext->FD[1] &= cl;	mynext->FD[0] &= cl;
				mynext->Guess4();
			}
		}
	}
}
void ZHONE::Guess6() {	//  solve digit 6/5
	if (zh1b_g.diag) {		
		cout << Char27out(zh1b_g.previous_ua_status[4]) << "guess 6" << endl;
		cout << Char27out(zh1b_g.upstream_unsolved_cells[4]) << "upstream unsolved" << endl;
		ImageCandidats();}
	if (!Update6()) return;
	int v = FD[5], vr = v >> 27;
	if (!vr) {
		uint32_t w6 = FD[5] & zh1b_g.fdsw[2][5];
		if((!zh1b_g.modegua)&& (!w6)) return;
		zh1b_g.upstream_unsolved_cells[3] =
			zh1b_g.upstream_unsolved_cells[4] & ~FD[5];
		if (zh1b_g.diag) {
			cout << Char27out(w6) << " guess6 ua exit solved" << endl;
		}
		zh1b_g.previous_ua_status[3] = w6 | zh1b_g.previous_ua_status[4];
		//				mynext->FD[4] &= (mynext->cells_unsolved | 7 << 27);
		uint32_t cl = (zh1b_g.upstream_unsolved_cells[3] | 7 << 27);
		FD[4] &= cl; FD[3] &= cl; FD[2] &= cl; FD[1] &= cl; FD[0] &= cl;
		Guess5();
		return;
	}
	vr = (-vr)&vr; // catch last bit
	v &= TblRowUnsolved[vr];// unknown in last unknown row
	if (zh1b_g.diag) cout << Char27out(v) << " row digit'6' for bit 0" << oct << vr << dec << endl;
	uint32_t cell;
	while (bitscanforward(cell, v)) {// first cell
		v ^= 1 << cell; // clear bit
		if (zh1b_g.diag) cout << " guess cell " << cell << endl;
		ZHONE * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Assign(5, cell); 
		if (mynext->UpdateDigit(5)) {// update digit 4/3
			if (zh1b_g.diag) {
				cout << "status after mynext->UpdateDigit(5)" << endl;
				cout << Char27out(mynext->FD[5]) << " FD[5] digit  "<< zh1b_g.digmap[5] + 1 << endl;
				mynext->ImageCandidats();
			}
			if (mynext->FD[5] >> 27) {
				mynext->Guess6();
			}
			else {// all digit 6 solved go to guess5 if not a zero deviation
				if (zh1b_g.diag)
					cout << Char27out(mynext->FD[5]) << " guess6 " << endl;
				uint32_t w6 = mynext->FD[5] & zh1b_g.fdsw[2][5];
				if (!w6) continue;
				zh1b_g.upstream_unsolved_cells[3] =
					zh1b_g.upstream_unsolved_cells[4] & ~mynext->FD[5];
				zh1b_g.previous_ua_status[3] = w6 | zh1b_g.previous_ua_status[4];
				uint32_t cl= (zh1b_g.upstream_unsolved_cells[3] | 7 << 27);
				mynext->FD[4] &= cl;	mynext->FD[3] &= cl;
				mynext->FD[2] &= cl;	mynext->FD[1] &= cl;	mynext->FD[0] &= cl;
				mynext->Guess5();
			}
		}
	}
}
void ZHONE::Guess7() {	//  solve digit 7/6
	if (zh1b_g.diag) {		cout << "guess 7" << endl;		ImageCandidats();	}
	if (!Update7()) return;;
	int v = FD[6], vr = v >> 27;
	if (!vr) {
		uint32_t w7 = FD[6] & zh1b_g.fdsw[2][6];
		if ((!zh1b_g.modegua) && (!w7)) return;
		zh1b_g.upstream_unsolved_cells[4] =
			zh1b_g.upstream_unsolved_cells[5] &~ FD[6];
		if (zh1b_g.diag)			cout << Char27out(w7) << " guess7 ua" << endl;
		zh1b_g.previous_ua_status[4] = w7 | zh1b_g.previous_ua_status[5];
		//				mynext->FD[4] &= (mynext->cells_unsolved | 7 << 27);
		uint32_t cl = (zh1b_g.upstream_unsolved_cells[4] | 7 << 27);
		FD[5] &= cl; FD[4] &= cl; FD[3] &= cl; 
		FD[2] &= cl; FD[1] &= cl; FD[0] &= cl;
		Guess6();
		return;
	}
	vr = (-vr)&vr; // catch last bit
	v &= TblRowUnsolved[vr];// unknown in last unknown row
	uint32_t cell;
	while (bitscanforward(cell, v)) {// first cell
		v ^= 1 << cell; // clear bit
		ZHONE * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Assign(6, cell);
		if (mynext->UpdateDigit(6)) {// update digit 4/3
			if (mynext->FD[6] >> 27) {
				mynext->Guess7();
			}
			else {// all digit 7 solved go to guess6 if not a zero deviation
				uint32_t w7 = mynext->FD[6] & zh1b_g.fdsw[2][6];
				if (!w7) continue;
				zh1b_g.upstream_unsolved_cells[4] =
					zh1b_g.upstream_unsolved_cells[5] & ~ mynext->FD[6];
				zh1b_g.previous_ua_status[4] = w7 | zh1b_g.previous_ua_status[5];
				uint32_t cl = (zh1b_g.upstream_unsolved_cells[4] | 7 << 27);
				mynext->FD[5] &= cl;	mynext->FD[4] &= cl;	mynext->FD[3] &= cl;
				mynext->FD[2] &= cl;	mynext->FD[1] &= cl;	mynext->FD[0] &= cl;
				mynext->Guess6();
			}
		}
	}
}
void ZHONE::AddMissingUAs(int * tcells, int ncells) {
	zh1b_g.type = 1;
	zh1b_g.ndigits = 9;
	zh1b_g.go_back = 0;
	int Digit_cell_Assigned[9];
	memset(Digit_cell_Assigned, 0, sizeof Digit_cell_Assigned);
	memset(this, 0, sizeof zhone[0]);
	InitOne_std_band();
	for (int icell = 0; icell < ncells; icell++) {
		int cell = tcells[icell],
			digit = zh1b_g.band0[cell],
			bit = 1 << cell;
		if ((FD[digit] & bit) == 0)  return;// check not valid entry
		Assign(digit, cell);
		Digit_cell_Assigned[digit] |= bit;
	}
	for (int i = 0; i < 9; i++)  FD[i] &=
		cells_unsolved |
		Digit_cell_Assigned[i] |
		07000000000;// don't touch rows unsolved
	//Debug(1);
	//ImageCandidats();
	if (FullUpdate()) Guess();
}

