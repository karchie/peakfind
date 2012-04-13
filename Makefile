all: 
	$(MAKE) -C src $@

clean:
	$(MAKE) -C src $@

tags:
	find src -name '*.[ch]' -print | xargs etags
