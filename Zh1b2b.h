
/*
ZhouSolver.h
Based on code posted to <http://forum.enjoysudoku.com/3-77us-solver-2-8g-cpu-testcase-17sodoku-t30470.html>
by user zhouyundong_2012.
The copyright is not specified.
*/

/* the BIT_SET_30 pattern is defined in a band for a digit
   it is a map of candidates plus 3 bits for unsolved rows
   bits 0-8 9-17 18-26 for the rows 27-29 for unsolved rows
*/
//#include "t_128GP.h"
// tables specific to the brute force located in zh4_tables
const extern int TblRowMask[8];// rows where single  found  000 to 111
const extern int  Tblstartblock[27]; // (i%3) * 27 in zhou brute force
const extern int TblShrinkMask[512];// existing minirows 000 to 111
const extern int TblComplexMask[512]; // keep mini rows still valid optimised process
const extern int TblMaskSingle[512]; // kill in other blocks locked column /box
const extern int TblMaskDouble[512];// kill for locked in box / column
const extern int TblColumnSingle[512]; // single in column applied to shrinked bloc
const extern int TblShrinkSingle[512]; // keep only rows with single
const extern int TblRowUniq[512]; // 1 is row not defined in block  mode  to 111
const extern T128 AssignMask_Digit[81];
extern uint64_t zh2b_start[20];

//const extern T128 AssignMask_OtherDigits[81];
// tables for UA collectors
extern int floors_2d[36];
extern int floors_3d[84]; // 84 3d  
extern int floors_4d[126];
extern int subsets_2d_3d[84][3];
extern int subsets_3d_4d[126][4];
extern int subsets_2d_4d[126][6];
extern int perm_0_8_size3[84][3];
extern int perm_0_8_size4[126][4];
//int * floors_3d_4d = floors_2d; // 84 3d   126 4d/5d
struct ZHOU;
/* class encapsulating the brute force with one known band
remaining clues are given in a 0-53 "cells" space
the 2 bands are located in a 64 bits field
unknown rows per digit are using 2x32 bits
first 5x6 bits for digits 0-4
second 4x6bits for digits 5-8

*/
struct ZH2B_GLOBAL { // global variables for the game table
	uint64_t * digsols; // pointer to solution grid per digit
	uint64_t ua_ret;
	BF64 val_init1_81;// , pairs, triplets;
	BF64 Digit_cell_Assigned_init[9];
	BF64 fd_sols[2][9];//start puzzle/ solution
	BF64 fd_revised[9];// gangster revision of the solution
	BF64 fdsw[3][9];//morphed digits puzzle/ solution rev
	BF64  myfd[9], cells_unsolved;// final Init5 status
	BF64 previous_ua_status[6];// index 0 is for digit 3
	uint64_t *tuaold, tua;// old floors uas /  fresh uas
	// ==============bands sols handling
	// sols collection ZH2B_1D
	BF64 *tsolw, *tuaw; // belong to the caller 
	BF64 mysol, mystart, andsol, myandsol;
	BF64 sols_buffer[3000], ua_buffer[3000];

	uint32_t nuaold, nua, ndigits;
	int tsd[7], ntsd, tsd2[7], ntsd2,
		socket_digits, isd1;// socket more
	int nsol, lim, icount, ntsol, single_applied, new_single_in_Update,
		rdigit, nctlg, go_back;
	// band UA collection active band pointers and UA table to build
	int modeguess;
	int  puz0[54], gangster[9];
	int nsolw;
	char * zsol, out54[55];// bands 12 in output mode
	char zdebug[82];

	ZH2B_GLOBAL();

	inline void InitsGetSols(int i,int ibuf){
		mystart = fdsw[1][i];
		mysol = fdsw[0][i];
		tsolw =&sols_buffer[ibuf];
		tuaw = &ua_buffer[ibuf];
	}
	void GetBands(int * g1, int * g2);
	void InitGangster(int * g0, int * g1);//common ot both GUA2s GUA3s
	uint64_t BuildUaret(BF64 * wsol);
};
/* 2 BF 128 per digit
	equivalent to F and Fcomp in Zhou
	Last 32 bits in FD[0] is a  bits field for unknown rows
	Last 32 bits in FD[1]] contain  digit mapped
*/
// class encapsulating the brute force 
struct ZH2B {// size 32 bytes 
	BF64 FD[9], CompFD[9], cells_unsolved, rows_unsolved;

	void Init_std_bands(); // after getbands in zh2b_g
	void Init_gang(); // after getbands in zh2b_g
	void DebugSol();
	inline void Assign(int digit, int cell, int xcell) {
		FD[digit] &= AssignMask_Digit[cell].u64[0];
		cells_unsolved.Clear(xcell);
		int ddig = 6 * digit;
		if (digit > 4) ddig += 2;// second bloc of 32 bits
		rows_unsolved.Clear(ddig + C_row[cell]);//6*digit + row
	}
	inline int Unsolved_Count() { return rows_unsolved.Count(); }
	
	void InitTclues(uint32_t * tclues, int n);
	void Init_2digits_banda(BF64  cellsbf);
	//void EndInit_2digits_bandb(int fl, int ibandb);

	uint64_t ValidXY(uint32_t * tclues, int n,int test=0);
	void DebugValidXY(uint32_t * tclues, int n, int test = 0);
	int Update();
	int FullUpdate();
	void ComputeNext();
	void GuessValidB12();
	void GuessGo(int dig, BF64 & wsol);

	int ApplySingleOrEmptyCells();
	char * SetKnown(char * zs);
	int Seta(int digit, int xcell);
	inline int SetaC(int digit, int cell) { return Seta(digit, C_To128[cell]); }
	inline void ClearCandidate(int dig, int xcell) { FD[dig].Clear(xcell); }
	int GetAllDigits(int cell);
	void Debug(int all=0);
	void ImageCandidats();
};

/* ZH2B5 same as ZH2B, limited to 5 digits
dedicated to UAs generati
*/
struct ZH2B5_GLOBAL { // global variables for the game table
	BF64 fdsw[3][5];//morphed digits puzzle/ solution rev
	BF64  myfd[5], cells_unsolved;// final Init5 status
	BF64 tuaf5[300];// collect uas for the cycle
	uint32_t nuaf5; // nuas found in the cycle
	int sizef5,// size limit for new uas
		single_applied, // loop control in full update
	    modevalid;// 0 base 1 gua mode
	uint32_t  ndigits,diag;
	//_______________________
	uint64_t FindUAsInit(int fl, int source = 1);
	void CollectUas5();//FindInit done
	void ValidSol5(uint64_t * sol);
};

struct ZH2B5 {
	BF64 FD[5], CompFD[5], cells_unsolved, rows_unsolved;
	inline void Assign(int digit, int cell, int xcell) {
		FD[digit] &= AssignMask_Digit[cell].u64[0];
		cells_unsolved.Clear(xcell);
		int ddig = 6 * digit;
		rows_unsolved.Clear(ddig + C_row[cell]);//6*digit + row
	}
	inline int Unsolved_Count() { return rows_unsolved.Count(); }
	//char * SetKnown(char * zs);
	int GetAllDigits(int cell);
	void ImageCandidats();

	//================== UA collector GUA mode
	int Update5();
	int FullUpdate5();
	int ApplySingleOrEmptyCells5();
	int Seta5(int digit, int xcell);
	inline void ComputeNext5() { if (FullUpdate5())Guess5(); }
	void Guess5();



};

/*ZH2B_1D class to solve one digit 2 bands
all valid solutions except the "all true" are stored in table
the table is supplied by the caller 
*/
struct ZH2B_1D_GLOBAL {
	BF64 *tsolw, *tuaw; // belong to the caller 
	BF64 mysol,  myandsol;
	int nsolw;
	//___________________
	int Go(BF64 & sol, BF64 & fde, BF64 *tsol, BF64 *tua,int ru);
};
struct ZH2B_1D {// size 32 bytes
	BF64 FD, CompFD;
	int GetSols(int ru);
	int GetAllSols(BF64 & fde, int ru, BF64 & fdsol);
	inline void Assign(int cell) {
		FD &= AssignMask_Digit[cell].u64[0];
	}
	void ComputeNext(int ru);
	int Update(int &ru);
	void Guess(int ru);
	int IsValid(uint64_t v);
};

 /* class encapsulating the brute force for one band
initial is a solved band giving the clues in each column
each bit map per digit is a 32 bit field
as in the original code, unknown rows per digit are located in bits 28-29-30 of the bit map
remaining clues are given in a 0-26 "cells" space

*/
struct ZHONE_GLOBAL { // global variables for the game table
	int diag;
	int nsol, lim, icount, ntsol, single_applied,
		new_single_in_Update,
		type,//gocatmode,    
		rdigit, nctlg, go_back;
	int pairs, triplets;
	char * zsol, out27[28];// band1 in output mode 

	// band UA collection active band pointers and UA table to build
	uint32_t *tua, nua;//   maximum 81  
	int * band0,*gangster;
	int floors_mini_row, digmap[9],digmap2[9];// to adjust pm
	uint32_t fd_sols[2][9];//start puzzle/ solution
	//uint32_t gua_gang_6_7[9]; // band initial pm for guas 6_7
	uint32_t fdsw[3][9];//morphed digits puzzle/ solution rev
	uint32_t ndigits,modegua;
	uint32_t previous_ua_status[6],// index 0 is for digit 3
		upstream_unsolved_cells[6];

	ZHONE_GLOBAL();
	void GetBand(int * b, uint32_t * t);
	void GetBand(uint32_t fd_solsb[2][9], int * b, uint32_t * t);
	inline void SetPat(char * pat, char * zsol, char * puzfinal) {
		pat = pat; zsol = zsol; puzfinal = puzfinal;
	};
	inline void InitIsvalid() { // usually after init 2 steps
		nsol = go_back = type = 0; lim = 1;
	}
	void ValidPuzzle(uint32_t * sol);
	void AddUA(uint32_t ua);
	void PrintTua();
	void FindMissingUAs();
};


 struct ZHONE {
	 uint32_t FD[9], CompFD[9], cells_unsolved;

	 // standard calls and basic process
	 inline void Assign(int digit, int cell) {
		 FD[digit] &= AssignMask_Digit[cell].u32[0];
		 cells_unsolved &= ~(1 << cell);
		 FD[digit] &= ~(1 << (27 + C_row[cell]));//9*digit + row
	 }
	 inline int Unsolved_Count() { return _popcnt32(cells_unsolved); }
	 int ApplySingleOrEmptyCells();
	 void InitOne_std_band(); // after getband in zh1b_g
	 void CheckSolPerDigit();     
	 int InitSudokux(GINT * t, int n);
	 void AddMissingUAs(int * tcells,int ncells);
	 int Update();
	 int Update4();
	 int Update6();
	 int Update7();
	 int FullUpdate();
	 void Guess();
	 inline void ComputeNext() { if (FullUpdate())	Guess(); }
	 char * SetKnown(char * zs);
	 void Seta(int digit, int xcell);
	 int Isvalid();
	 int GetAllDigits(int cell);
	 void Debug(int all = 0);
	 void DebugDigit(int digit);
	 void ImageCandidats();
	 //==================================================
	 // uacollector 
	 void Checkstart();
	 int  Start_nFloors(int floors);
	 void Start_Uas_Mini(int floors, int floors_mini_row);
	 void ApplyGangsterChanges(int * g0, int * g1);
	 void InitGuess();
	 int UpdateDigit(int digit);
	 void Guess2();	void Guess3();	void Guess4();
	 void Guess5(); void Guess6();  void Guess7();
	 void Set2(int cell);
 };

