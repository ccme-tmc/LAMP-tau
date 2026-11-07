all: manual.pdf

manual.pdf: manual.tex
	xelatex -interaction=nonstopmode manual.tex
	xelatex -interaction=nonstopmode manual.tex

clean:
	rm -f manual.aux manual.log manual.out manual.toc
