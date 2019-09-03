# sk_chk19
search puzzles 17-19 clues in a solution grid

finding low clues puzzles in a solution grid is done for long in gridchecker.
This is an alternative process derived from the search of 17 clues puzzles through a scan of all solution grids.
As in sk_s17, another repository from me, the process is derived from blue's proposal for the search of 17s (blue, an active member of the "New Sudoku Player's Forum")
Here, having three bands A,B,C 
a) All valid bands A of size <= x are scanned
b) for each valid band A all bands B giving a valid {bands A+B) are produced
c) then valid bands C giving a global valid puzzle (unique and minimal) are produced.

Note, the project starts September 2019. The first repository is only valid for the UAs generation. The rest of the code is just copied from sk_s17. A first version of the program is expected in October 2019

