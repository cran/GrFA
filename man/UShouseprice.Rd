\name{UShouseprice}
\alias{UShouseprice}
\docType{data}
\title{
Housing price data for 16 states in the U.S.
}
\description{
This dataset contains the Zillow Home Value Index (ZHVI) at the county level for single-family residences and condos with 1, 2, 3, 4, or 5+ bedrooms. It focuses on the middle tier of home values (33rd to 67th percentile) and features smoothed, seasonally adjusted values presented on a monthly basis. The data spans 16 U.S. states from January 2000 to April 2023. Within each state, the data is organized as a matrix, and the data for all states is compiled into a list.
}
\usage{data("UShouseprice")}
\format{
The dataset is structured as a list containing 16 elements, with each element corresponding to a state. Each element is a matrix where the columns represent time series data for house prices at the county level. Each time series has a length of 280, representing monthly data points from January 2000 to April 2023. The number of columns in each matrix varies, ranging from 90 to 250, depending on the number of counties and bedroom categories in the state. The columns are labeled with the county name and bedroom count (e.g., “Pulaski County bd1” for one-bedroom homes or “Garland County bd5” for homes with five or more bedrooms). This structure provides a comprehensive and organized representation of the Zillow Home Value Index (ZHVI) across multiple counties and bedroom categories for the 16 states included in the dataset.
}
\source{
The original data is downloaded from the website of Zillow.
}
\details{
The column names of the data matrix represent county names combined with bedroom counts. For example, "Pulaski County bd1" indicates the house price in Pulaski County for one-bedroom homes, while "Garland County bd5" refers to the house price in Garland County for homes with more than five bedrooms.

The abbreviations and full names of these 16 states are as follows:

AR: Arkansas

CA: California

CO: Colorado

FL: Florida

GA: Georgia

KY: Kentucky

MD: Maryland

MI: Michigan

NC: North Carolina

NJ: New Jersey

NY: New York

OH: Ohio

OK: Oklahoma

PA: Pennsylvania

TN: Tennessee

VA: Virginia
}
\examples{
data(UShouseprice)
log_diff = function(x){
  T = nrow(x)
  res = log(x[2:T,]/x[1:(T-1),])*100
  scale(res, center = TRUE, scale = TRUE)
}
UShouseprice1 = lapply(UShouseprice, log_diff)
}
\keyword{datasets}
