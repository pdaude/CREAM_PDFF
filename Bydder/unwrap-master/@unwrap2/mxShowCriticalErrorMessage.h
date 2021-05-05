// Beautified mex error message macro.
//
// Inspired by
//   http://www.advanpix.com/2016/02/14/short-and-informative-error-messages-from-mex/
//

// Use macro to expand __FILE__ and __LINE__ correctly
#define mxShowCriticalErrorMessage(...) err_fn(basename(__FILE__),__LINE__,##__VA_ARGS__)

// Use overloads to handle __VA_ARGS__ correctly
void err_fn(const char *fn_name, int line_no, const char *err_message, int err_code)
{
    const int nargs = 5;
    mxArray *err_args[nargs];
    err_args[0] = mxCreateString("\n%s(%i): %s (%i).\n");
    err_args[1] = mxCreateString(fn_name);
    err_args[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    err_args[3] = mxCreateString(err_message);
    err_args[4] = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(err_args[2]) = line_no;
    *mxGetPr(err_args[4]) = err_code;
    mexCallMATLAB(0,0,nargs,err_args,"error");
}

void err_fn(const char *fn_name, int line_no, const char *err_message)
{
    const int nargs = 4;
    mxArray *err_args[nargs];
    err_args[0] = mxCreateString("\n%s(%i): %s.\n");
    err_args[1] = mxCreateString(fn_name);
    err_args[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    err_args[3] = mxCreateString(err_message);
    *mxGetPr(err_args[2]) = line_no;
    mexCallMATLAB(0,0,nargs,err_args,"error");
}

void err_fn(const char *fn_name, int line_no, int err_code)
{
    err_fn(fn_name, line_no, "Error occurred", err_code);
}

void err_fn(const char *fn_name, int line_no)
{
    err_fn(fn_name, line_no, "Error occurred");
}
