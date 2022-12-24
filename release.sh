clang main.c glad_gl.c -I inc -Ofast -lglfw -lm -o wiggle
strip --strip-unneeded wiggle
upx --lzma --best wiggle