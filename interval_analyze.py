import pandas as pd
import plotly.express as px


def create_intervals(haplotype_dict):
    """
    This function will create a list for each child in the following format:
    [(interval number is the index) {start position: , end position: ,
    haplotype: (1 or 2)]
    each interval starts with the position of a variant from haplotype 1 or 2,
    and ends when the next variant is from the opposite haplotype, where a new
    interval will start
    """
    intervals = []
    current_interval = None

    for position, value in haplotype_dict.items():
        chromosome_num = int(haplotype_dict[position][-3])
        cur_haplotype = value[-2]
        if current_interval is None:
            # Start a new interval
            current_interval = {"start": position, "end": position,
                                "haplotype": cur_haplotype,
                                "chromosome": chromosome_num}
        elif cur_haplotype == current_interval["haplotype"]:
            # Continue the current interval
            current_interval["end"] = position
            # Adding the chromosome of the interval
            current_interval["chromosome"] = chromosome_num
        else:
            # Start a new interval as haplotype changed
            intervals.append({"start": current_interval["start"],
                              "end": current_interval["end"],
                              "haplotype": current_interval["haplotype"],
                              "chromosome": chromosome_num})
            current_interval = {"start": position, "end": position,
                                "haplotype": cur_haplotype}

    # Add the last interval
    if current_interval is not None:
        intervals.append({"start": current_interval["start"],
                          "end": current_interval["end"],
                          "haplotype": current_interval["haplotype"],
                          "chromosome": current_interval["chromosome"]})

    return intervals


def shared_interval(interval_lists):
    """
    This function will create a new list containing intervals that are shared
    in all the lists given, according to the haplotype
    """
    # Initialize shared_intervals with the intervals from the first list
    shared_intervals = interval_lists[0]

    # Iterate through the remaining lists
    for interval_list in interval_lists[1:]:
        # A temporary list to store shared intervals for the current list
        temp_shared_intervals = []

        # Iterate through each interval in the current list
        for interval_1 in shared_intervals:
            for interval_2 in interval_list:
                check = interval_1["chromosome"]
                if (
                        interval_1["haplotype"] == interval_2["haplotype"]
                        and interval_1["start"] <= interval_2["end"]
                        and interval_1["end"] >= interval_2["start"]
                        and interval_1["chromosome"] == interval_2["chromosome"]
                ):
                    # Calculate the intersection of intervals
                    start = max(interval_1["start"], interval_2["start"])
                    end = min(interval_1["end"], interval_2["end"])
                    temp_shared_intervals.append({"start": start, "end": end,
                                                  "haplotype": interval_1[
                                                      "haplotype"],
                                                  "chromosome": interval_1[
                                                      "chromosome"]})

        # Update shared_intervals with the current shared intervals
        shared_intervals = temp_shared_intervals

    return shared_intervals


def plot_interval(interval_list, plot_title):
    """
    This function plots intervals as straight lines, where each interval is
    represented by a line starting from the "start" to "end" keys on the
    x-axis,
    and the height determined by the "haplotype" key on the y-axis.
    """
    # Create a DataFrame for Plotly Express
    data = {"Start": [], "End": [], "Haplotype": [], "Interval": []}

    for i, interval in enumerate(interval_list):
        start_position = interval["start"]
        end_position = interval["end"]
        haplotype = interval["haplotype"]

        # Append data to the DataFrame
        data["Start"].extend([start_position, end_position])
        data["End"].extend([start_position, end_position])
        data["Haplotype"].extend([haplotype, haplotype])
        data["Interval"].extend([f'Interval {i + 1}', f'Interval {i + 1}'])

    # Create a DataFrame
    df = pd.DataFrame(data)

    # Create an interactive line plot using Plotly Express
    fig = px.line(df, x="Start", y="Haplotype", color="Interval",
                  labels={'Start': 'Chromosome Position',
                          'Haplotype': 'Haplotype'},
                  title=plot_title)

    # Show the plot in an HTML window
    fig.show()
