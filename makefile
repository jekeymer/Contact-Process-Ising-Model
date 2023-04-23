all:
	gcc "CPIM v2.0.c" mt64.c -lm -o "CPIM v2.0" `pkg-config --cflags gtk+-3.0` `pkg-config --libs gtk+-3.0`
