# Local Genome Browser — third-party runtime

The local genome browser uses:

- **IGV.js 3.8.5** — MIT license. It is not bundled in this ZIP; the browser downloads
  the pinned `igv.min.js` once from the official jsDelivr npm package and caches it in
  `reference_data/browser_runtime/`.
  Source: `https://cdn.jsdelivr.net/npm/igv@3.8.5/dist/igv.min.js`
  License: `https://github.com/igvteam/igv.js/blob/master/LICENSE`
  Citation: Robinson JT, Thorvaldsdottir H, Turner D, Mesirov JP.
  *igv.js: an embeddable JavaScript implementation of the Integrative Genomics Viewer*.
  Bioinformatics. 2023;39(1):btac830.

- **Flask 3.1.3** — installed only when menu item 24 is launched if it is absent.
  The application binds only to `127.0.0.1`.

No genomic file is uploaded by the browser. BAM/BAI, VCF/TBI and FASTA/FAI are served
from the local pipeline directory to the local browser through an allowlisted localhost endpoint.
