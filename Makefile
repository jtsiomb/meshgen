PREFIX = /usr/local

src = $(wildcard src/*.cc)
obj = $(src:.cc=.o)
dep = $(obj:.o=.d)
inst_hdr = src/mesh.h src/meshgen.h src/geom.h src/object.h
bin = meshgen

CXXFLAGS = -pedantic -Wall -g -Isrc -DPREFIX=\"$(PREFIX)\" -fPIC
LDFLAGS = -lGL -lGLU -lglut -lGLEW -lm -lgmath -ldl -rdynamic

$(bin): $(obj)
	$(CXX) -o $@ $(obj) $(LDFLAGS)

-include $(dep)

%.d: %.cc
	@$(CPP) $(CXXFLAGS) $< -MM -MT $(@:.d=.o) >$@

.PHONY: clean
clean:
	rm -f $(obj) $(bin)

.PHONY: cleandep
cleandep:
	rm -f $(dep)

.PHONY: install
install: $(bin)
	mkdir -p $(DESTDIR)$(PREFIX)/bin
	mkdir -p $(DESTDIR)$(PREFIX)/include/meshgen
	cp $(bin) $(DESTDIR)$(PREFIX)/bin
	cp $(inst_hdr) $(DESTDIR)$(PREFIX)/include/meshgen

.PHONY: uninstall
uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/$(bin)
	rm -f $(DESTDIR)$(PREFIX)/include/meshgen/*.h
	rmdir $(DESTDIR)$(PREFIX)/include/meshgen
