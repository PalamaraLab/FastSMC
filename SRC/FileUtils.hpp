#ifndef FILEUTILS_HPP
#define FILEUTILS_HPP

#include <vector>
#include <string>
#include <fstream>

#include "StringUtils.hpp"

#include <boost/iostreams/filtering_stream.hpp>

namespace FileUtils {

void openOrExit(std::ifstream &stream, const std::string &file,
                std::ios_base::openmode mode = std::ios::in);

void openWritingOrExit(std::ofstream &stream, const std::string &file,
                       std::ios_base::openmode mode = std::ios::out);

void requireEmptyOrReadable(const std::string &file);

void requireEachEmptyOrReadable(const std::vector <std::string> &fileList);

void requireEmptyOrWriteable(const std::string &file);

std::vector <std::string> parseHeader(const std::string &fileName,
                                      const std::string &delimiters);

int lookupColumnInd(const std::string &fileName, const std::string &delimiters,
                    const std::string &columnName);

double readDoubleNanInf(std::istream &stream);

std::vector < std::pair <std::string, std::string> > readFidIids(const std::string &file);

class AutoGzIfstream {
  boost::iostreams::filtering_istream boost_in;
  std::ifstream fin;
public:
  static int lineCount(const std::string &file);

  void openOrExit(const std::string &file, std::ios_base::openmode mode = std::ios::in);
  void close();
  template <class T> AutoGzIfstream& operator >> (T &x) {
    boost_in >> x;
    return *this;
  }
  operator bool() const;
  AutoGzIfstream& read(char *s, std::streamsize n);
  int get();
  double readDoubleNanInf();
  void clear();
  AutoGzIfstream& seekg(std::streamoff off, std::ios_base::seekdir way);
  friend AutoGzIfstream& getline(AutoGzIfstream& in, std::string &s);
};

AutoGzIfstream& getline(AutoGzIfstream& in, std::string &s);

class AutoGzOfstream {
  boost::iostreams::filtering_ostream boost_out;
  std::ofstream fout;
public:
  void openOrExit(const std::string &file, std::ios_base::openmode mode = std::ios::out);
  void close();
  template <class T> AutoGzOfstream& operator << (const T &x) {
    boost_out << x;
    return *this;
  }
  AutoGzOfstream& operator << (std::ostream & (*manip)(std::ostream&));
  void unsetf(std::ios_base::fmtflags);
  operator bool() const;
};

}

#endif
