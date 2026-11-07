# files = gauss_elimination linsys_solver LU_matrix LU_solver matrix_print QR_matrix QR_solver
linsys_solver_files = linsys_solver.o gauss_elimination.o matrix_print.o LU_matrix.o LU_solver.o QR_matrix.o QR_solver.o gauss_seidel.o



linsys_solver: $(linsys_solver_files)
	gcc -g -Wall $(linsys_solver_files) -o linsys_solver -lm


linsys_solver.o: linsys_solver.c
	gcc -c -g -Wall $<
gauss_elimination.o: gauss_elimination.c
	gcc -c -g -Wall $<
LU_matrix.o: LU_matrix.c
	gcc -c -g -Wall $<
LU_solver.o: LU_solver.c
	gcc -c -g -Wall $<
matrix_print.o: matrix_print.c
	gcc -c -g -Wall $<
QR_matrix.o: QR_matrix.c
	gcc -c -g -Wall $<
QR_solver.o: QR_solver.c
	gcc -c -g -Wall $<
gauss_seidel.o: gauss_seidel.c
	gcc -c -g -Wall $<

clean:
	rm -f *.o