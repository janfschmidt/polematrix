TEXCMD = pdflatex --shell-escape
BIBCMD = biber

all: polematrix.pdf

polematrix.pdf: polematrix.tex polem_defs.sty polem.bib coord.pdf
	$(TEXCMD) polematrix
	$(BIBCMD) polematrix
	$(TEXCMD) polematrix

clean:
	rm -rf auto/
	rm -rf _minted-polematrix/
	rm -f *.aux *.bbl *.bcf *.blg *.log *.out *.run.xml *.toc
	rm polematrix.pdf
