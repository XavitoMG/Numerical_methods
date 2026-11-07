
zero_solver_files = finite_difference.o eval_derivative.o vector_field.o newton_raphson.o LU_matrix.o LU_solver.o zero_solver.o

zero_solver: $(zero_solver_files)
	gcc -g -Wall $(zero_solver_files) -o zero_solver -lm

finite_difference.o: finite_difference.c
	gcc -c -g -Wall $<
eval_derivative.o: eval_derivative.c
	gcc -c -g -Wall $<
vector_field.o: vector_field.c
	gcc -c -g -Wall $<
newton_raphson.o: newton_raphson.c
	gcc -c -g -Wall $<
LU_matrix.o: LU_matrix.c
	gcc -c -g -Wall $<
LU_solver.o: LU_solver.c
	gcc -c -g -Wall $<
zero_solver.o: zero_solver.c
	gcc -c -g -Wall $<


clean:
	rm -f *.o	
