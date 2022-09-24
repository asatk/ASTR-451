CC = gcc
CFLAGS = -Wall -Wextra -Werror -g
LDFLAGS = -lm

main: parse.c main.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

clean:
	$(RM) main *.o *~ core.[1-9]*