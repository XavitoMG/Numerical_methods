BIN = bin
SRC_IVP = src/initial_value_problem
SRC_LS = src/linear_systems
SRC_NLO = src/non_linear_optimization
SRC_BVP = src/boundary_value_problem
APP = apps
INC = include

ivp_solver_files =                 \
    src/vector_field.o             \
	$(SRC_IVP)/rk78.o              \
	$(SRC_IVP)/flow.o

ivp_solver_main = $(APP)/ivp_solver.o

linsys_solver_files =              \
	$(SRC_LS)/gauss_elimination.o  \
    $(SRC_LS)/gauss_seidel.o       \
	$(SRC_LS)/LU_matrix.o          \
	$(SRC_LS)/LU_solver.o          \
	$(SRC_LS)/matrix_print.o       \
	$(SRC_LS)/QR_matrix.o          \
	$(SRC_LS)/QR_solver.o

linsys_solver_main = $(APP)/linsys_solver.o

zero_solver_files =                \
	src/vector_field.o             \
	$(SRC_NLO)/finite_difference.o \
    $(SRC_NLO)/eval_derivative.o   \
	$(SRC_NLO)/newton_raphson.o    \
	$(SRC_LS)/matrix_print.o       \
	$(SRC_LS)/LU_matrix.o          \
	$(SRC_LS)/LU_solver.o          \

zero_solver_main = $(APP)/zero_solver.o

bvp_solver_files =                 \
	$(zero_solver_files)           \
	$(ivp_solver_files)            \
	$(SRC_BVP)/shooting_method.o   \
	$(SRC_BVP)/shooting_function.o

bvp_solver_main = $(APP)/bvp_solver.o

	

all: $(BIN)/ivp_solver $(BIN)/linsys_solver $(BIN)/zero_solver $(BIN)/bvp_solver

$(BIN)/ivp_solver: $(ivp_solver_files) $(ivp_solver_main)
	gcc -g -Wall $^ -o $@ -lm
$(BIN)/linsys_solver: $(linsys_solver_files) $(linsys_solver_main)
	gcc -g -Wall $^ -o $@ -lm
$(BIN)/zero_solver: $(zero_solver_files)  $(zero_solver_main)
	gcc -g -Wall $^ -o $@ -lm
$(BIN)/bvp_solver: $(bvp_solver_files) $(bvp_solver_main)
	gcc -g -Wall $^ -o $@ -lm

$(SRC_IVP)/%.o: $(SRC_IVP)/%.c $(INC)/*.h
	gcc -c -g -Wall $< -o $@ -I$(INC)
$(SRC_LS)/%.o: $(SRC_LS)/%.c $(INC)/*.h
	gcc -c -g -Wall $< -o $@ -I$(INC)
$(SRC_NLO)/%.o: $(SRC_NLO)/%.c $(INC)/*.h
	gcc -c -g -Wall $< -o $@ -I$(INC)
$(SRC_BVP)/%.o: $(SRC_BVP)/%.c $(INC)/*.h
	gcc -c -g -Wall $< -o $@ -I$(INC)
$(APP)/%.o: $(APP)/%.c $(INC)/*.h
	gcc -c -g -Wall $< -o $@ -I$(INC)
src/vector_field.o: src/vector_field.c
	gcc -c -g -Wall $< -o $@ -I$(INC)

clean:
	rm -f $(SRC_IVP)/*.o $(SRC_LS)/*.o $(SRC_NLO)/*.o $(SRC_BVP)/*.o $(APP)/*.o src/*.o $(BIN)/*
