#' Housing price data for 16 states in the U.S.
#'
#' @description
#' This dataset contains the Zillow Home Value Index (ZHVI) at the county level for single-family residences and condos with 1, 2, 3, 4, or 5+ bedrooms. It focuses on the middle tier of home values (33rd to 67th percentile) and features smoothed, seasonally adjusted values presented on a monthly basis. The data spans 16 U.S. states from January 2000 to April 2023. Within each state, the data is organized as a matrix, and the data for all states is compiled into a list.
#'
#' @format The dataset is structured as a list containing 16 elements, with each element corresponding to a state. Each element is a matrix where the columns represent time series data for house prices at the county level.
#' \describe{
#'   \item{Rows}{Each time series has a length of 280, representing monthly data points from January 2000 to April 2023.}
#'   \item{Columns}{The number of columns in each matrix varies, ranging from 90 to 250, depending on the number of counties and bedroom categories in the state.}
#'   \item{Labels}{The columns are labeled with the county name and bedroom count (e.g., "Pulaski County bd1" for one-bedroom homes or "Garland County bd5" for homes with five or more bedrooms).}
#' }
#'
#' @details
#' The column names of the data matrix represent county names combined with bedroom counts. For example, "Pulaski County bd1" indicates the house price in Pulaski County for one-bedroom homes, while "Garland County bd5" refers to the house price in Garland County for homes with more than five bedrooms.
#'
#' The abbreviations and full names of these 16 states are as follows:
#' \itemize{
#'   \item AR: Arkansas
#'   \item CA: California
#'   \item CO: Colorado
#'   \item FL: Florida
#'   \item GA: Georgia
#'   \item KY: Kentucky
#'   \item MD: Maryland
#'   \item MI: Michigan
#'   \item NC: North Carolina
#'   \item NJ: New Jersey
#'   \item NY: New York
#'   \item OH: Ohio
#'   \item OK: Oklahoma
#'   \item PA: Pennsylvania
#'   \item TN: Tennessee
#'   \item VA: Virginia
#' }
#'
#' @source The original data is downloaded from the website of Zillow.
#'
#' @references
#' Aggregated Projection Method: A New Approach for Group Factor Model. Jiaqi Hu, Ting Li, Xueqin Wang (2025). Journal of the American Statistical Association, doi:10.1080/01621459.2025.2491154
#'
#' @usage data(UShouseprice)
#'
#' @examples
#' data(UShouseprice)
#'
#' # Helper function to calculate log differences and scale
#' log_diff <- function(x) {
#'   T <- nrow(x)
#'   res <- log(x[2:T, ] / x[1:(T - 1), ]) * 100
#'   scale(res, center = TRUE, scale = TRUE)
#' }
#'
#' # Apply to all states
#' UShouseprice1 <- lapply(UShouseprice, log_diff)
#'
#' @docType data
#' @keywords datasets
"UShouseprice"
