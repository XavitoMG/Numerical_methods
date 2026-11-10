BIN = bin
SRC_IVP = src/initial_value_problem
SRC_LS = src/linear_systems
SRC_NLO = src/non_linear_optimization
APP = apps
INC = include

ivp_solver_files =                 \
    src/vector_field.o             \
	$(SRC_IVP)/rk78.o              \
	$(APP)/ivp_solver.o

linsys_solver_files =              \
	$(SRC_LS)/gauss_elimination.o  \
    $(SRC_LS)/gauss_seidel.o       \
	$(SRC_LS)/LU_matrix.o          \
	$(SRC_LS)/LU_solver.o          \
	$(SRC_LS)/matrix_print.o       \
	$(SRC_LS)/QR_matrix.o          \
	$(SRC_LS)/QR_solver.o          \
	$(APP)/linsys_solver.o

zero_solver_files =                \
	src/vector_field.o             \
	$(SRC_NLO)/finite_difference.o \
    $(SRC_NLO)/eval_derivative.o   \
	$(SRC_NLO)/newton_raphson.o    \
	$(SRC_LS)/LU_matrix.o          \
	$(SRC_LS)/LU_solver.o          \
	$(APP)/zero_solver.o

all: $(BIN)/ivp_solver $(BIN)/linsys_solver $(BIN)/zero_solver

$(BIN)/ivp_solver: $(ivp_solver_files)
	gcc -g -Wall $^ -o $@ -lm
$(BIN)/linsys_solver: $(linsys_solver_files)
	gcc -g -Wall $^ -o $@ -lm
$(BIN)/zero_solver: $(zero_solver_files)
	gcc -g -Wall $^ -o $@ -lm

$(SRC_IVP)/%.o: $(SRC_IVP)/%.c $(INC)/*.h
	gcc -c -g -Wall $< -o $@ -I$(INC)
$(SRC_LS)/%.o: $(SRC_LS)/%.c $(INC)/*.h
	gcc -c -g -Wall $< -o $@ -I$(INC)
$(SRC_NLO)/%.o: $(SRC_NLO)/%.c $(INC)/*.h
	gcc -c -g -Wall $< -o $@ -I$(INC)
$(APP)/%.o: $(APP)/%.c $(INC)/*.h
	gcc -c -g -Wall $< -o $@ -I$(INC)
src/vector_field.o: src/vector_field.c
	gcc -c -g -Wall $< -o $@ -I$(INC)

clean:
	rm -f $(SRC_IVP)/*.o $(SRC_LS)/*.o $(SRC_NLO)/*.o $(APP)/*.o src/*.o $(BIN)/*
