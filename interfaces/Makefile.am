AUTOMAKE_OPTIONS = foreign

clean-local::
	rm -f interfaces.h

BUILT_SOURCES = interfaces.h

interfaces.h:interfaces.h.in Makefile
	python $(abs_top_srcdir)/pre-process-defs.py $@ $<
