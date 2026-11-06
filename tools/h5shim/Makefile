libh5_deflate_shim.so: h5_deflate_shim.c
	$(CC) -O2 -Wall -shared -fPIC -o $@ $< -ldl
clean:
	rm -f libh5_deflate_shim.so
