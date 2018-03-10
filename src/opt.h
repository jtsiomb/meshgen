#ifndef OPT_H_
#define OPT_H_

struct Options {
	char *fname;
	char *cxx;
};
extern Options opt;

bool parse_options(int argc, char **argv);

#endif	// OPT_H_
