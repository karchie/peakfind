all: lib

lib:
	$(MAKE) -C src $@

tags:
	find src -name '*.[ch]' -print | xargs etags
