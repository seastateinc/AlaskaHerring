EXEC=ham
TPLS=../../src/$(EXEC).tpl
DATA=$(EXEC).dat
PATH=../../src/build
ARGS=

ifdef DEBUG
	DIST=$(PATH)/debug/$(EXEC)
else
	DIST=$(PATH)/release/$(EXEC)
endif


all: $(EXEC) $(EXEC).rep

$(EXEC): $(DIST)
	cp $(DIST) $@



$(EXEC).rep: $(EXEC) 
	./$(EXEC) -ind $(DATA) $(ARGS)



run: $(DATA)
	./$(EXEC) -ind $(DATA) $(ARGS)


clean:
	@rm -f $(EXEC) 
	@rm -f admodel.cov
	@rm -f admodel.hes
	@rm -f admodel.dep
	@rm -f fmin.log
