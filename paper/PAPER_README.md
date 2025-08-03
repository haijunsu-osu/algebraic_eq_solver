# Scientific Paper: A Robust Numerical Solver for Trigonometric Algebraic Systems

This directory contains a comprehensive scientific paper describing the trigonometric equation solver.

## Paper Status
- ✅ LaTeX source completed (`paper.tex`)
- ⏳ PDF compilation pending (requires LaTeX installation)
- 📋 Ready for arXiv submission

## Directory Contents
- **`paper.tex`**: Complete LaTeX source with arXiv template
- **`compile_paper.bat`**: Windows compilation script
- **`compile_paper.sh`**: Unix/Mac/Linux compilation script
- **`COMPILE_INSTRUCTIONS.md`**: Detailed setup and compilation guide
- **`PAPER_README.md`**: This documentation file

## Quick Compilation

### Windows
```cmd
cd paper
compile_paper.bat
```

### Unix/Mac/Linux
```bash
cd paper
chmod +x compile_paper.sh
./compile_paper.sh
```

### Online (No LaTeX installation needed)
1. Upload `paper.tex` to [Overleaf](https://www.overleaf.com/)
2. Compile online and download PDF
3. Copy the PDF to the root directory as `paper.pdf`

## Paper Structure
1. **Title & Author**: Haijun Su (haijunsu-osu)
2. **Abstract**: Comprehensive overview of the solver capabilities
3. **Introduction**: Problem motivation and literature context
4. **Method**: Complete mathematical formulation including:
   - Generic quartic polynomial approach
   - Singular matrix handling (zero and rank-1 cases)
   - Unified algorithm description
5. **Numerical Examples**: 
   - Non-singular case with two solutions
   - Zero matrix case (free parameter)
   - Rank-1 case with four solutions
   - Performance analysis
6. **Conclusions**: Robustness and practical applications
7. **References**: 9 relevant academic sources

## Workflow
1. **Compile**: Run compilation script or use online LaTeX
2. **Verify**: Check that `paper.pdf` is generated
3. **Commit**: The PDF will be copied to root directory for git tracking
4. **Submit**: Upload to arXiv with GitHub repository link

## Next Steps
1. Compile `paper.tex` to generate `paper.pdf`
2. Commit the PDF to repository: `git add paper.pdf && git commit -m "Add compiled scientific paper PDF"`
3. Submit to arXiv with GitHub repository link
4. Add arXiv link to repository README

The paper provides a complete scientific description suitable for academic publication and serves as comprehensive documentation for the solver implementation.
