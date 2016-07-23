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
	@echo "All done"


$(EXEC).rep: $(EXEC) 
	./$(EXEC) -ind $(DATA) $(ARGS)
	@echo "Done here to"

run: $(DATA)
	./$(EXEC) -ind $(DATA) $(ARGS)


clean:
	rm -f $(EXEC) 
