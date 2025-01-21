# Variables
CC = gcc
CFLAGS = -Wall -Wextra -Wshadow -Wconversion -pedantic -std=c11 -O2 -fstack-protector-strong -D_FORTIFY_SOURCE=2
LIB_NAME = libcalc.a
SRC = src/add.c src/subtract.c
OBJ = src/add.o src/subtract.o

.PHONY: all clean

all: $(LIB_NAME)

$(LIB_NAME): $(OBJ)
	ar rcs $@ $^

src/%.o: src/%.c include/%.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(LIB_NAME)
