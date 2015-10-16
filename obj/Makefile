# Edmanuel Torres 21/10/2004

include ../makeinclude

all: $(PRGS)

$(PRGS):
	$(CXX) -o $(XDIR)$@ $(PRGS).o $(OBJS) $(LFLAGS)

clean:
	rm -rf *o $(PRGS)