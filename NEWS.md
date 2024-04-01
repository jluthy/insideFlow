# insideFlow 1.1.0

*This version was released on 31Mar2024.*

## Minor Changes

- Changed the default behavior of `loadFlowJoCSV` to automatically detect the `cellId` and `eventNumber` columns if they are present in the input data. This change simplifies the function call for most users, as the `cellId` and `eventNumber` parameters are now optional. Changed to use `fread` from `data.table` for faster CSV reading.

- Added additional parameters `cellId` and `eventNumber` to `loadFlowJoCSV` for flexible specification of cell ID and event number columns.

- Improved color palette customization in `plotFlowSOMResults`, including support for `RColorBrewer` and `viridis` color palettes.

- Added tests for `loadFlowJoCSV` function and creation of insideFlow objects.

# insideFlow 1.0.1

*This version was released on 28May2023.*

## Bug Fixes

- Fixed an issue where missing `cellId` or `eventNumber` columns would halt execution. Now, missing columns are handled gracefully, with warnings issued and the relevant data structures left empty.

# insideFlow 1.0.0

*This version was released on 27May2023.*

- Initial release of the insideFlow package, providing foundational tools and methods for analyzing FlowJo workspace data within R.
- Features include data import from FlowJo, basic data manipulation and preprocessing, and integration with FlowSOM for advanced clustering and visualization.
