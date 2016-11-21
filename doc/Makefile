prj=avrescale

ver=0.0

figgps=$(wildcard fig/*.gp)
pdffigs=$(patsubst %.gp,%.pdf, $(figgps))
epsfigs=$(patsubst %.gp,%.eps, $(figgps))

texpdf = $(prj).pdf $(prj)_notes.pdf

alltexpdf = $(patsubst %.tex,%.pdf, $(wildcard $(prj)*.tex))

regular: $(texpdf)

all: $(alltexpdf)

$(alltexpdf) : %.pdf : %.tex simul.bib $(pdffigs) Makefile
	make -C fig
	pdflatex $<
	bibtex $*
	pdflatex $<
	pdflatex $<

$(pdffigs): $(epsfigs) $(figgps)
	make -C fig

clean:
	$(RM) -rf *.aux *.bbl *.blg *.dvi *.log *Notes.bib *~
	rstrip.py -Rlv

zip: $(prj)docV$(ver).zip

# make a zip for the current version

$(prj)docV$(ver).zip::
	mkdir -p old/V$(ver)
	rsync -avz $(prj).tex old/V$(ver)/$(prj)_annotatedV$(ver).tex
	rsync -avz Makefile old/V$(ver)/
	rsync -avz README.md old/V$(ver)/
	rsync -avz deannote.py old/V$(ver)/
	rsync -avz rmbiburl.py old/V$(ver)/
	rsync -avz checkeq.py old/V$(ver)/
	python rmbiburl.py -i simul.bib -o old/V$(ver)/simul.bib
	cd old/V$(ver) && \
	  ln -sf ../../fig .
	cd old/V$(ver) && \
	  python deannote.py -v -c -i $(prj)_annotatedV$(ver).tex -o $(prj)V$(ver).tex
	make all -C old/V$(ver)/
	cd old/V$(ver) && \
	  zip -r ../../$(prj)docV$(ver).zip *.tex *.pdf simul.bib *.py \
	       $(epsfigs) $(pdffigs) $(figgps) Makefile README.md

#	rsync -avz $(prj)_notes.tex old/V$(ver)/$(prj)_notesV$(ver).tex


$(prj)doc.zip: $(prj).tex $(prj).pdf \
  	       $(prj)_notes.tex $(prj)_notes.pdf \
  	       simul.bib \
	       $(epsfigs) $(pdffigs) $(figgps) \
	       Makefile README.md
	zip -r $@ $^

Bossman: zip
	cp $(prj)docV$(ver).zip /Bossman/cz1/avrescale/doc/
	#rsync -avL $(prj)docV$(ver).zip cz1@129.109.88.204:/Bossman/cz1/avrescale/doc/
