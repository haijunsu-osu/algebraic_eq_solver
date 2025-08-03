@echo off
REM Compile the scientific paper and prepare for repository inclusion

echo Compiling scientific paper...

REM Check if pdflatex is available
where pdflatex >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo Error: pdflatex not found. Please install a LaTeX distribution:
    echo   - Windows: MiKTeX or TeX Live
    echo   - Online: Use Overleaf at https://www.overleaf.com/
    pause
    exit /b 1
)

REM Compile the paper (run twice for cross-references)
echo Running pdflatex (first pass)...
pdflatex paper.tex

echo Running pdflatex (second pass)...
pdflatex paper.tex

REM Clean up intermediate files
echo Cleaning up intermediate files...
if exist paper.aux del paper.aux
if exist paper.log del paper.log
if exist paper.out del paper.out
if exist paper.fls del paper.fls
if exist paper.fdb_latexmk del paper.fdb_latexmk
if exist paper.synctex.gz del paper.synctex.gz

REM Check if PDF was created successfully
if exist paper.pdf (
    echo ✅ Paper compiled successfully: paper.pdf
    echo 📋 Ready to commit PDF to repository
    echo 🚀 Ready for arXiv submission
    
    REM Show file info
    dir paper.pdf
) else (
    echo ❌ Compilation failed. Check the LaTeX source for errors.
    pause
    exit /b 1
)

pause
