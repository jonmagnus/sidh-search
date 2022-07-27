
PDFS = $(foreach tuple,\
       	4_3 5_4 6_5, \
	$(foreach agent,alice bob,\
	$(agent)_$(tuple).pdf \
))

all : $(PDFS)

%.pdf : %.dot
	sfdp -Tpdf -Goverlap=scale $< -o $@

define GEN_RULE
$(1)_$(2)_$(3).dot :
	../explore_graph $$@ $2 $3 $(if $(subst $1,,alice), 0, 1)
endef

$(foreach src, $(PDFS),\
	$(eval $(call GEN_RULE,$(word 1, $(subst _, ,$(src))),$(word 2, $(subst _, ,$(src))),$(word 3, $(subst ., ,$(subst _, ,$(src)))))))

clean:
	rm *.dot *.pdf
