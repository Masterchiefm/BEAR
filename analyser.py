import markdown
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import numpy as np
import pandas as pd
from io import BytesIO
import base64
import tabulate
from Bio import pairwise2


class SangerBaseCall:
    def __init__(self, sanger_file):
        self.sanger_data = SeqIO.read(sanger_file, "abi")
        self.sanger_seq = str(self.sanger_data.seq).upper()

        self.labels = {}
        signal_length = self.sanger_data.annotations["abif_raw"]["PLOC1"][-1]
        base_locations = self.sanger_data.annotations["abif_raw"]["PLOC1"]
        self.annotation = {}
        self.location_base_annotation = {}

        # 列出每个碱基位置上分别是什么碱基和他对应的测序峰值位置。
        for i in range(len(self.sanger_seq)):
            try:
                self.annotation[i] = [self.sanger_seq[i], base_locations[i]]
                self.location_base_annotation[base_locations[i]] = self.sanger_seq[i]
            except:
                continue

        # 列出每个信号列中的信号值
        channels = {"G": "DATA9", "A": "DATA10", "T": "DATA11", "C": "DATA12"}  # GATC
        self.trace = {}
        for c in channels:
            self.trace[c] = self.sanger_data.annotations["abif_raw"][channels[c]]

        # 根据信号峰值填入碱基位置
        for i in range(signal_length):
            if i in base_locations:
                self.labels[i] = self.location_base_annotation[i]
            else:
                self.labels[i] = ""

        # 整体的噪音
        self.noise = {"G": [], "A": [], "T": [], "C": []}

    def locTartet(self, target_seq, max_mismatch=10):
        # Find the match site and location
        target_is_reversed = False
        target_seq = str(target_seq).upper()
        alignment_1 = pairwise2.align.localms(target_seq, self.sanger_seq, 2, -.1, -4, -2, one_alignment_only=True)
        target_seq_reversed = str(Seq(target_seq).reverse_complement())
        alignment_2 = pairwise2.align.localms(target_seq_reversed, self.sanger_seq, 2, -.1, -4, -2,
                                              one_alignment_only=True)
        try:
            alignment_1[0]
        except:
            alignment_1 = [[0, 0, 0]]
        try:
            alignment_2[0]
        except:
            alignment_2 = [[0, 0, 0]]

        if alignment_2[0][2] > alignment_1[0][2]:
            alignment = alignment_2
        else:
            alignment = alignment_1
        # print(alignment)

        if target_is_reversed:
            # print('reverse')
            location_start = alignment[0][3]
            location_end = alignment[0][4]
            # print(location_start,location_end)
            match_seq = str(alignment[location_end][location_start]).upper()
        else:
            location_start = alignment[0][3]
            location_end = alignment[0][4]
            match_seq = str(alignment[0][1][location_start:location_end]).upper()

        return (match_seq, location_start, location_end)

    def plotTarget(self, target_seq, plot_window=10):
        from matplotlib import pyplot as plt
        plt.cla()
        plt.close("all")
        # plt.figure(dpi=150)
        # plt.figure(figsize=(20, 4))
        plt.rcParams['figure.figsize'] = (15, 4)
        colors = ("black", "green", "red", "blue")  # GATC
        target_location = self.locTartet(str(target_seq))
        print(target_location)

        if target_location[1] < target_location[2]:
            target_reverse = False
            base_start = target_location[1] - plot_window
            base_end = target_location[2] + plot_window
        else:
            target_reverse = True
            base_start = target_location[2] - plot_window
            base_end = target_location[1] + plot_window

        if base_start * base_end <= 0:
            print("Plot window size out of range")

        plot_start = self.annotation[base_start][1]
        plot_end = self.annotation[base_end][1]

        # GATC trace
        # self.trace["G"]
        trace_G = self.trace["G"][plot_start:plot_end]
        trace_A = self.trace["A"][plot_start:plot_end]
        trace_T = self.trace["T"][plot_start:plot_end]
        trace_C = self.trace["C"][plot_start:plot_end]

        labels = list(self.labels.values())[plot_start: plot_end]

        fig, ax = plt.subplots()

        for base_i in range(len(labels)):  # 在图上标注碱基
            base = labels[base_i]
            if base != "":
                plt.text(base_i - 1, -100, base)

        plt.tick_params(axis="x",
                        length=0,  # tick长度
                        width=0,  # tick的宽度
                        labelsize=0,  # 刻度值大小
                        bottom=False

                        )
        # plt.xticks(ticks=range(plot_end - plot_start), labels=labels, label="sequence")
        plt.xlabel("Sequence")
        plt.ylabel('Sanger Signal Intensity')
        G = ax.plot(trace_G, color=colors[0], label='G')
        A = ax.plot(trace_A, color=colors[1], label='A')
        T = ax.plot(trace_T, color=colors[2], label='T')
        C = ax.plot(trace_C, color=colors[3], label='C')

        sg_plot_start = self.annotation[target_location[1]][1]
        sg_plot_end = self.annotation[target_location[2]][1]

        if sg_plot_start < sg_plot_end:
            sg_plot_end = self.annotation[target_location[2] - 1][1]
            SG = ax.plot([sg_plot_start - plot_start, sg_plot_end - plot_start], [-150, -150], color="brown",
                         label='Target', marker=">")



        else:
            sg_plot_start = self.annotation[target_location[1]][1]
            SG = ax.plot([sg_plot_end - plot_start, sg_plot_start - plot_start - 9], [-150, -150], color="purple",
                         label='Target', marker="<")

        plt.legend()
        return plt

    def getTargetData(self, target_seq, annalyse_window=10):
        """获取靶位点上各个碱基的信号占比，以及靶位点外的噪音平均值"""
        # plt.rcParams['figure.figsize'] = (15, 4)
        # colors = ("black", "green", "red", "blue")  # GATC
        target_location = self.locTartet(str(target_seq))
        # print(target_location)

        if str(target_location) == "None":
            return "None"

        if target_location[1] < target_location[2]:
            target_reverse = False
            base_start = target_location[1] - annalyse_window
            base_end = target_location[2] + annalyse_window
            base_5p_loc = target_location[1]
            # print(base_start, base_end)
        else:
            target_reverse = True
            base_start = target_location[2] - annalyse_window
            base_end = target_location[1] + annalyse_window
            base_5p_loc = target_location[1]
            # print(base_start, base_end)

        if base_start * base_end <= 0:
            print("Plot window size out of range")

        plot_start = self.annotation[base_start][1]
        plot_end = self.annotation[base_end][1]

        # 噪声统计区间
        noise_start_1 = 100
        noise_end_1 = base_start - 20
        noise_start_2 = base_end + 20

        try:
            noise_end_2 = 700
        except Exception as e:
            print("sanger file is too short!")
            noise_end_2 = base_end + 120

        self.noise = {"G": [], "A": [], "T": [], "C": []}  # 用来记录噪音，后续计算平均值和方差
        # print(noise_start_1, noise_end_1)
        for i in range(noise_start_1, noise_end_1):
            if i in self.location_base_annotation:
                peak = i
                peak_base = self.location_base_annotation[peak]
                peak_base_num = list(self.location_base_annotation.keys()).index(peak)
                peak_left_border = peak - 3
                peak_right_border = peak + 3

                signal = {}
                signal["G"] = self.trace["G"][peak_left_border:peak_right_border]
                signal["A"] = self.trace["A"][peak_left_border:peak_right_border]
                signal["T"] = self.trace["T"][peak_left_border:peak_right_border]
                signal["C"] = self.trace["C"][peak_left_border:peak_right_border]

                for base in ["G", "A", "T", "C"]:
                    # print(base+peak_base)
                    if base != peak_base:
                        # print(signal[base])
                        self.noise[base].append(signal[base])

        for i in range(noise_start_2, noise_end_2):
            if i in self.location_base_annotation:
                peak = i
                peak_base = self.location_base_annotation[peak]
                peak_base_num = list(self.location_base_annotation.keys()).index(peak)
                peak_left_border = peak - 3
                peak_right_border = peak + 3

                signal = {}
                signal["G"] = self.trace["G"][peak_left_border:peak_right_border]
                signal["A"] = self.trace["A"][peak_left_border:peak_right_border]
                signal["T"] = self.trace["T"][peak_left_border:peak_right_border]
                signal["C"] = self.trace["C"][peak_left_border:peak_right_border]

                for base in ["G", "A", "T", "C"]:
                    # print(base+peak_base)
                    if base != peak_base:
                        # print(signal[base])
                        self.noise[base].append(signal[base])

        result = {}
        # fig, ax = plt.subplots()
        sheet = pd.DataFrame(index=["G", "A", "T", "C", "sequence", "guide"],dtype=object)
        for i in range(plot_start, plot_end):
            if i in self.location_base_annotation:
                peak = i
                peak_base = self.location_base_annotation[peak]
                peak_base_num = list(self.location_base_annotation.keys()).index(peak)
                peak_trace = self.trace[peak_base]
                # print(peak_trace)eak_val = peak_trace[peak]

                """
                注释掉这一段。这段代码用于查找整个信号的汉宁窗口，然后再以该窗口计算信号、噪音。
                但是，窗口边界可能与其它碱基的窗口重叠，从而导致噪音误判。

                只采用信号峰值附近的应该就够了。
                peak_left_data = np.array(peak_trace[peak-10:peak][::-1])
                peak_right_data = np.array(peak_trace[peak:peak+10])


                peak_left_distance = argrelextrema(peak_left_data,comparator=np.less_equal, order=3)
                peak_right_distance = argrelextrema(peak_right_data,comparator=np.less_equal, order=3)
                peak_left_width = (peak_left_distance[0][0])
                peak_right_width = (peak_right_distance[0][0])
                peak_left_border = peak - peak_left_width
                peak_right_border = peak + peak_right_width

                """

                peak_left_border = peak - 3
                peak_right_border = peak + 3

                signal_G = self.trace["G"][peak_left_border:peak_right_border]
                signal_A = self.trace["A"][peak_left_border:peak_right_border]
                signal_T = self.trace["T"][peak_left_border:peak_right_border]
                signal_C = self.trace["C"][peak_left_border:peak_right_border]

                signal_intensity = {}
                signal_intensity["G"] = sum(signal_G)
                signal_intensity["A"] = sum(signal_A)
                signal_intensity["T"] = sum(signal_T)
                signal_intensity["C"] = sum(signal_C)
                signal_intensity["ALL"] = sum(
                    [signal_intensity["G"], signal_intensity["A"], signal_intensity["T"], signal_intensity["C"]])

                if target_reverse:
                    # print("base_5p_loc",base_5p_loc)
                    peak_base_relative_num = base_5p_loc - peak_base_num

                    if peak_base_relative_num <= 0:
                        peak_base_relative_num = base_5p_loc - peak_base_num - 1
                else:
                    peak_base_relative_num = peak_base_num - base_5p_loc
                    if peak_base_relative_num >= 0:
                        peak_base_relative_num = peak_base_relative_num + 1

                result[peak_base_relative_num] = [peak_base, peak_base_num, signal_intensity]

        # plot_sequence = pd.DataFrame()
        for location in result:
            base = result[location][0]
            sheet.loc["G", location] = (result[location][2]["G"] / result[location][2]["ALL"] * 100)
            sheet.loc["A", location] = (result[location][2]["A"] / result[location][2]["ALL"] * 100)
            sheet.loc["T", location] = (result[location][2]["T"] / result[location][2]["ALL"] * 100)
            sheet.loc["C", location] = (result[location][2]["C"] / result[location][2]["ALL"] * 100)
            sheet.loc["sequence", location] = base

            if (location < len(target_seq)) and (location > 0) and (target_reverse == False):
                sheet.loc["guide", location] = target_seq[location - 1]
            elif (location < len(target_seq)) and (location > 0) and (target_reverse == True):
                sheet.loc["guide", location] = target_seq[location - 1]

            else:
                sheet.loc["guide", location] = ""

        return sheet

    def plotTargetHeatmap(self, data_sheet, base_of_interest="", color_map="Blues"):
        from matplotlib import pyplot as plt
        plt.cla()
        plt.close("all")
        sequence = data_sheet.loc["sequence"]
        columns = ["G", "A", "T", "C", "Guide"]
        sg_data = []

        # 给sgRNA的位置也赋值，从而带上颜色
        for base in data_sheet.loc["guide"]:
            if base == "":
                sg_data.append(-15)
            else:
                sg_data.append(-50)

        data = np.array([data_sheet.loc["G"],
                         data_sheet.loc["A"],
                         data_sheet.loc["T"],
                         data_sheet.loc["C"],#sg_data
                         ]).astype("float")

        # fig, axes = plt.subplots(2,1,figsize=(20, 6))
        fig, ax = plt.subplots(figsize=(20, 6))
        # ax = axes[0]
        ax.set_title(
            "Signal intensity of each base in the sequence " + "".join(list(data_sheet.loc["guide"])))
        im = ax.imshow(data, cmap=color_map)

        ax.set_xticks(np.arange(len(sequence)), labels=sequence)
        ax.set_yticks(np.arange(len(columns)-1), labels=columns[:-1])

        ax.tick_params(top=False, bottom=False, labeltop=False, labelbottom=True, left=False)
        # ax2 = axes[1]
        # ax2.imshow([sg_data])
        # ax2.set_xticks(np.arange(len(sequence)), labels=sequence)


        for i in range(len(columns)):
            for j in range(len(sequence)):
                if i < 4:
                    text = ax.text(j, i, round(data[i, j], 1), ha="center", va="center", color="black")
                else:
                    text = ax.text(j, i + 0.2, list(sequence.index)[j], ha="center", va="center", color="red")
        if base_of_interest == "":
            pass
        elif base_of_interest == "G":
            y = 0
            ax.add_patch(plt.Rectangle((-0.5, y - 0.5), len(sequence), 1, fill=False, edgecolor="red", lw=5))
        elif base_of_interest == "A":
            y = 1
            ax.add_patch(plt.Rectangle((-0.5, y - 0.5), len(sequence), 1, fill=False, edgecolor="red", lw=5))
        elif base_of_interest == "T":
            y = 2
            ax.add_patch(plt.Rectangle((-0.5, y - 0.5), len(sequence), 1, fill=False, edgecolor="red", lw=5))
        elif base_of_interest == "C":
            y = 3
            ax.add_patch(plt.Rectangle((-0.5, y - 0.5), len(sequence), 1, fill=False, edgecolor="red", lw=5))

        return plt

    def plt2md(self, plt):
        buffer = BytesIO()
        plt.savefig(buffer, format='png')

        # 转换base64并以utf8格式输出
        pic_base64 = base64.b64encode(buffer.getvalue()).decode('utf8')
        md = """![](data:image/png;base64,""" + pic_base64 + """)"""
        return md

    def generateReport(self, target_seq, report_name, save_path, plot_window=5, base_of_interest=""):
        sanger_plot = self.plotTarget(target_seq=target_seq, plot_window=10)
        sanger_plot_md = self.plt2md(sanger_plot)

        perc_around_guide_sheet = self.getTargetData(target_seq=target_seq, annalyse_window=plot_window)
        heatmap_plot = self.plotTargetHeatmap(perc_around_guide_sheet, base_of_interest)
        heatmap_plot_md = self.plt2md(heatmap_plot)

        title = "# BaseEdit Report of " + report_name + "\n"
        sep = "\n\n---------"
        description = """## Description 描述
        \n\nthis file contains the sanger sequencing report around the guide 
        """ + target_seq + f"\n\n 此文件包含了靶序列({target_seq})以及其上下游部分序列的测序分析结果。若有碱基编辑，则有套峰出现。"

        sub_title_1 = "\n\n## Sanger sequencing result around " + target_seq + "\n\n"
        Fig1_annotation = "\n\nThis figure shows the Sanger sequencing signal quality." + f"此图展示了靶位点及其附件的Sanger测序质量"

        sub_title_2 = "\n\n## Base call result in each position around " + target_seq + "靶位点及其上下游各碱基的信号强度\n"
        Fig2_annotation = """\n\nThis figure shows the possibilities of each base at each sequencing position. 
                                Guide position at the botton line shows the relative position to your guide RNA.\n\n此图展示了各位点上各个碱基的信号强度"""
        if base_of_interest != "":
            Fig2_annotation = Fig2_annotation + "\n\nThe red rectangle is the base you want to convert to.\n\n红框为你要改变成为的碱基。"
            description = description + "\n\nThis sample is supposed to be converted to " + base_of_interest + f"\n\n此样品期望被转为{base_of_interest}"

        sub_title_3 = "\n\n## Base call raw data in each position around " + target_seq + "原始数据表\n"
        sheet = "\n\n" + perc_around_guide_sheet.to_markdown() + "\n\n"
        sheet_annotation = """\n\n This sheet shows the possibilities of each base at each sequencing position.\n\nIt's the same as the heatmap above.
                                \n\nYou can copy the entire sheet and pasted it to Excel for further analysis.\n\n原始数据表,利用这些数据可以生成上图。如果需要进行更多的分析或者重新绘图，请直接复制下表并粘贴到excel进行后续处理\n\n
                                """
        content = [title, description, sep, sub_title_1, sanger_plot_md, Fig1_annotation, sep, sub_title_2,
                   heatmap_plot_md, Fig2_annotation, sep, sub_title_3, sheet_annotation, sheet]

        html_content = markdown.markdown("".join(content), extensions=["tables"])
        with open(save_path, "w") as f:
            header = """
<!DOCTYPE html>
<head>
	<meta charset="UTF-8">
</head>
<html>
	<body>"""
            footer = """</body>
</html>"""
            f.write(html_content)


if __name__ == "__main__":

    a = SangerBaseCall("example.ab1")
    a.generateReport("CACTGGAATGACACACGCCC", "report_name", "test.html",plot_window=5)
