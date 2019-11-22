#!/bin/bash

for FILE in _images/*.pdf; do
    ls "${FILE}"
    gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dNOPAUSE -dQUIET -dBATCH  -dAutoRotatePages=/None -dUseCIEColor -dCompressFonts=true -r150 -sOutputFile="${FILE}.1" "${FILE}"
    mv "${FILE}.1" "${FILE}"
done
