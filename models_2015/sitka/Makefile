EXEC=ham
TPLS=../../src/$(EXEC).tpl
DATA=$(EXEC).dat
SRCS=../../src/build
ARGS=

ifdef DEBUG
	DIST=$(SRCS)/debug/$(EXEC)
else
	DIST=$(SRCS)/release/$(EXEC)
endif


all: $(EXEC) $(EXEC).rep

$(EXEC): $(DIST)
	ln -s $(DIST) $@



$(EXEC).rep: $(EXEC) $(DATA)
	./$(EXEC) -ind $(DATA) $(ARGS)



run: $(DATA)
	./$(EXEC) -ind $(DATA) $(ARGS)


clean:
	@rm -f $(EXEC) 
	@rm -f admodel.cov
	@rm -f admodel.hes
	@rm -f admodel.dep
	@rm -f fmin.log
	@rm -f $(EXEC).rep $(EXEC).par ham.[rpls]
