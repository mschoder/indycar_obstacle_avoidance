#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>
/*
 * A class to read data from a csv file, customized for "trackfile_parsed.csv"
 */
class CSVReader
{
    std::string fileName;
    std::string delimeter;
public:
    CSVReader(std::string filename, std::string delm = ",") :
            fileName(filename), delimeter(delm)
    { }
    // Function to fetch data from a CSV File
    std::vector<std::vector<double> > getData();
};
/*
* Parses through csv file line by line and returns the data
* in vector of vector of strings.
*/
std::vector<std::vector<double> > CSVReader::getData()
{
    std::ifstream file(fileName);
    std::vector<std::vector<double> > dataList;
    std::string line = "";
    // Iterate through each line and split the content using delimeter
    while (getline(file, line))
    {
        if (line[0] == '#') {continue;} // skip first line in csv file
        std::vector<std::string> vec;
        boost::algorithm::split(vec, line, boost::is_any_of(delimeter));
        std::vector<double> dvec(vec.size());
        std::transform(vec.begin(), vec.end(), dvec.begin(), [](const std::string& val)
        {
            return std::stod(val);
        });
        dataList.push_back(dvec);
    }
    // Close the File
    file.close();
    return dataList;
}