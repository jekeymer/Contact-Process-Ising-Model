all:
	gcc CPIM.c mt64.c -lm -o CPIM `pkg-config --cflags gtk+-3.0` `pkg-config --libs gtk+-3.0`
