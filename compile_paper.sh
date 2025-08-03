#!/bin/bash
# Compile the scientific paper and prepare for repository inclusion

echo "Compiling scientific paper..."

# Check if pdflatex is available
if ! command -v pdflatex &> /dev/null; then
    echo "Error: pdflatex not found. Please install a LaTeX distribution:"
    echo "  - Windows: MiKTeX or TeX Live"
    echo "  - macOS: MacTeX"
    echo "  - Linux: texlive-full package"
    exit 1
fi

# Compile the paper (run twice for cross-references)
echo "Running pdflatex (first pass)..."
pdflatex paper.tex

echo "Running pdflatex (second pass)..."
pdflatex paper.tex

# Clean up intermediate files
echo "Cleaning up intermediate files..."
rm -f paper.aux paper.log paper.out paper.fls paper.fdb_latexmk paper.synctex.gz

# Check if PDF was created successfully
if [ -f "paper.pdf" ]; then
    echo "✅ Paper compiled successfully: paper.pdf"
    echo "📋 Ready to commit PDF to repository"
    echo "🚀 Ready for arXiv submission"
    
    # Show file size
    echo "PDF size: $(ls -lh paper.pdf | awk '{print $5}')"
else
    echo "❌ Compilation failed. Check the LaTeX source for errors."
    exit 1
fi
