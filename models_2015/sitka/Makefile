EXEC=ham
DAT=sitka_2.dat
CTL=sitka_2.ctl

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


# Retrospective Rule
NRET=1
retdirs := $(shell echo 'cat(formatC(1:$(NRET),digits=3,flag="0"))' | R --vanilla --slave)
makdirs := $(foreach dir, $(retdirs), $(dir)/makdirs)
rundirs := $(foreach dir, $(retdirs), $(dir)/rundirs)

$(makdirs): $(EXEC)
	mkdir $(@D)
	cp $(EXEC) $(DATA) $(DAT) $(CTL) Makefile $(@D)
	cd $(@D); touch makdirs

$(rundirs): 
	-cd $(@D); make run ARGS='-retro $(@D)'
	cd $(@D); touch rundirs

retro: $(EXEC) $(makdirs) $(rundirs)
	touch RETRODONE

# Simulation Rule

clean:
	@rm -f $(EXEC) 
	@rm -f admodel.cov
	@rm -f admodel.hes
	@rm -f admodel.dep
	@rm -f fmin.log
	@rm -f $(EXEC).rep $(EXEC).par ham.[rpls]
	@rm -r 0* 