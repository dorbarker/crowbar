import sys
import functools

def user_msg(*messages):
    """Wrapper for print() that prints to stderr"""
    print(*messages, file=sys.stderr)


def logtime(name):
    """Function decorator that print to stderr the runtime of the
    decorated function.
    """

    def decorator(func):
        """Interface between wrapper and the outer logtime() function"""

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            """Wraps func and prints to STDERR the runtime of func"""

            msg = 'Elapsed time for {}: {}'

            before = datetime.now()

            result = func(*args, *kwargs)

            after = datetime.now()

            user_msg(msg.format(name, after - before))

            return result
        return wrapper
    return decorator
