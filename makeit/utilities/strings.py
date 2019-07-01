
def string_or_range_to_float(text):
    """Translate a number or range in string format to a float.

    If the string represents a range, return the average of that range.

    Args:
        text (str): Number or range string to be converted to a float.

    Returns:
        float or None: Number converted / average of range or None if it could
            not be converted.
    """

    try:
        return float(text)
    except Exception as e:
        if text.count('-') == 1:  # 20 - 30
            try:
                x = text.split('-')
                return (float(x[0]) + float(x[1])) / 2.0
            except Exception as e:
                print(e)
        elif text.count('-') == 2:  # -20 - 0
            try:
                x = text.split('-')
                return (-float(x[0]) + float(x[1])) / 2.0
            except Exception as e:
                print(e)
        elif text.count('-') == 3:  # -20 - -10
            try:
                x = text.split('-')
                return (-float(x[0]) - float(x[1])) / 2.0
            except Exception as e:
                print(e)
        else:
            print(e)
    return None
