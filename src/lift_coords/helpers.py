def file_exists_nonempty(file_path):
    return True if os.path.isfile(file_path) \
                   and os.path.getsize(file_path) > 0 else False


def page(sliceable, step=500, start=0):
    """Yield slices of sliceable object.

    Args:
        step (int): number of items per slice.
        start (int): starting index for first slice.
    """
    nxt = start
    past_end = False
    while not past_end:
        last, nxt = nxt, nxt + step
        yield sliceable[last:nxt]
        if nxt > len(sliceable):
            past_end = True
