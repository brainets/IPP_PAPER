VPATH=.
OBJDIR=.obj
DEPDIR=.dep
CXX=c++
#
INCLUDES+=-I.
LDLIBS+=-lpthread -lpng
ifeq ($(DEBUG),)
DEBUG=no
endif
ifeq ($(DEBUG),no)
	CPUFLAGS=-O6 -fPIC -m64 -funroll-loops -std=c++11
else
	CPUFLAGS=-O0 -fPIC -m64 -ggdb -rdynamic -std=c++11
        #CPUFLAGS=-O1 -fPIC -m64 -ggdb -rdynamic -std=c++11
endif
# rules
$(OBJDIR)/%.o:%.cpp
	@if ! [ -d $(OBJDIR) ]; then mkdir $(OBJDIR); echo objects in $(OBJDIR)...; fi
	$(CXX) $(CXXFLAGS) -c $< -o $(OBJDIR)/$(@F)

$(DEPDIR)/%.d:%.cpp
	@if ! [ -d $(DEPDIR) ]; then mkdir $(DEPDIR); echo dependencies in $(DEPDIR)...; fi
	@set -e; rm -f $(DEPDIR)/$(@F);\
	$(CXX) $(CXXFLAGS) -MM -MT $(OBJDIR)/$(*F).o $< > $(DEPDIR)/$(@F).$$$$; \
	sed 's,\($(*F)\)\.o[ :]*,\1.o $(DEPDIR)/$(@F) : ,g' < $(DEPDIR)/$(@F).$$$$ > $(DEPDIR)/$(@F); \
	rm -f $(DEPDIR)/$(@F).$$$$
# end rules
WARNINGS=-Wall
CXXFLAGS=$(CPUFLAGS) $(WARNINGS) $(INCLUDES)
#
# objects
#
OBJ+=nmring.o util.o threadpool.o #hdf5util.o
OBJECTS=$(addprefix $(OBJDIR)/,$(OBJ))
DEPS   =$(addprefix $(DEPDIR)/,$(OBJ:.o=.d))
#
ONERINGOBJ=oneRing_model.o
ONERINGOBJECTS=$(addprefix $(OBJDIR)/,$(ONERINGOBJ))
ONERINGDEPS   =$(addprefix $(DEPDIR)/,$(ONERINGOBJ:.o=.d))

TWORINGOBJ=fb2rings_model.o
TWORINGOBJECTS=$(addprefix $(OBJDIR)/,$(TWORINGOBJ))
TWORINGDEPS   =$(addprefix $(DEPDIR)/,$(TWORINGOBJ:.o=.d))

THREERINGOBJ=ff3rings_model.o
THREERINGOBJECTS=$(addprefix $(OBJDIR)/,$(THREERINGOBJ))
THREERINGDEPS   =$(addprefix $(DEPDIR)/,$(THREERINGOBJ:.o=.d))

DEPS+=$(ONERINGDEPS) 
all: $(DEPS) oneRing_model fb2rings_model ff3rings_model

# dependencies
-include $(DEPS)

#
# targets
#
oneRing_model: $(OBJECTS) $(ONERINGOBJECTS)
	$(CXX) $(CXXFLAGS) $(LFLAGS) $^ -o $@ $(LDLIBS)

fb2rings_model: $(OBJECTS) $(TWORINGOBJECTS)
	$(CXX) $(CXXFLAGS) $(LFLAGS) $^ -o $@ $(LDLIBS)

ff3rings_model: $(OBJECTS) $(THREERINGOBJECTS)
	$(CXX) $(CXXFLAGS) $(LFLAGS) $^ -o $@ $(LDLIBS)	

# phony targets
tar:
	./autotar

clean:
	rm -rf $(DEPDIR) $(OBJDIR)
	rm -f oneRing_model fb2rings_model ff3rings_model ./*~

