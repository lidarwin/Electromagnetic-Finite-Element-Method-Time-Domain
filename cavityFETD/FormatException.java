package cavityFETD;

/** An unchecked exception that represents any kind of user error in the
 *  input to the formatter.
 *  @author Darwin Li
 */
class FormatException extends RuntimeException {

    /** A FormatException with no message. */
    FormatException() {
    }

    /** A FormatException for which .getMessage() is MSG. */
    FormatException(String msg) {
        super(msg);
    }

    /** Returns an exception containing an error message formatted according
     *  to FORMAT and ARGS, as for printf or String.format. Typically, one uses
     *  this by throwing the result in a context where there is a 'try' block
     *  that handles it by printing the message (esp. via reportError). */
    static FormatException error(String format, Object... args) {
        return new FormatException(String.format(format, args));
    }

    /** Print error message formed from arguments FORMAT and ARGS, whose
     *  meaning is as for printf. */
    static void reportError(String format, Object... args) {
        System.err.printf(format, args);
        System.err.println();
        _totalErrors += 1;
        System.exit(1);
    }

    /** Returns the total number of calls to reportError. */
    static int getTotalErrors() {
        return _totalErrors;
    }

    /** Cumulative errors encountered.  Assumes that 'error' is called to
     *  report each error. */
    private static int _totalErrors;

}
