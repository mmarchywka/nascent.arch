

static char *
is_good_dvi_file(const char *filename, Boolean from_history)
{
    static char canonical_path[MAXPATHLEN + 1];
    Boolean tried_dvi_extension = False;
    /* following allocates real_filename */
    char *real_filename = find_dvi_file(filename, &tried_dvi_extension, from_history);
    char *ret;
    FILE *f = NULL;
    dviErrFlagT errflag;
        
    if (real_filename == NULL)
        return NULL;

    if ((ret = REALPATH(real_filename, canonical_path)) == NULL) {
        /* REALPATH failed, use real_filename */
        strncpy(canonical_path, real_filename, MAXPATHLEN);
        canonical_path[MAXPATHLEN] = '\0';
        ret = canonical_path;
    }
    free(real_filename);

    /* check for correct DVI files */
    if ((f = XFOPEN(ret, OPEN_MODE)) != NULL) {
        TRACE_EVENTS((stderr, "watching: new file opened successfully."));
        if (process_preamble(f, &errflag)
            && find_postamble(f, &errflag) 
            && read_postamble(f, &errflag, False
#if DELAYED_MKTEXPK
                              , False
#endif
                              )) {
            fclose(f);
            return ret;
        }
        fclose(f);
        if (!from_history)
            XDVI_FATAL((stderr, "%s: %s.", filename, get_dvi_error(errflag)));
        return NULL;
    }
    else { 
        if (!from_history)
            XDVI_FATAL((stderr, "Could not open `%s': %s.", filename, strerror(errno)));
        return NULL;
    }
}

