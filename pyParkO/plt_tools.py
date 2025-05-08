from abc import ABCMeta, abstractmethod
from matplotlib import pyplot as plt
import numpy as np
import os

plt.rcParams['font.sans-serif'] = ['Times New Roman']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams.update({'font.size': 12})
file_path = os.path.dirname(__file__)


# Abstract class for plot tools
class PlotC(metaclass=ABCMeta):
    @abstractmethod
    def plt_conterf(self, matrix, x_ticks=None, y_ticks=None, title=None, colorbar=None, name=None, padding_dict={}):
        """
        plot filled contour map direct from matrix
        :param matrix: n × p
        :param x_ticks: list, e.g. [-10, -5, 0, 5, 10]
        :param y_ticks: list, e.g. [-10, -5, 0, 5, 10]
        :param title: str if None, no title
        :param colorbar:  None, vertical or horizontal
        :param name: str if None, plt.show else plt.save
        :param padding_dict: left < right  bottom < top        float value adjust the position of the graph
            :parameter left: float default 0.125
            :parameter bottom: float default 0.1
            :parameter right: float default 0.9
            :parameter top: float default 0.9
        :return: None
        """
        pass

    # 3D map
    @abstractmethod
    def plt_3d(self, matrix, x_ticks=None, y_ticks=None, title=None, name=None, padding_dict={}):
        """
        plot 3d-map
        :param mtrix: n * p
        :param x_ticks: list, e.g. [-10, -5, 0, 5, 10]
        :param y_ticks: list, e.g. [-10, -5, 0, 5, 10]
        :param title: str if None, no title
        :param name: if None, plt.show else plt.save
        :param padding_dict: left < right  bottom < top        float value adjust the position of the graph
            :parameter left: float default 0.125
            :parameter bottom: float default 0.1
            :parameter right: float default 0.9
            :parameter top: float default 0.9
        :return: None
        """
        pass

    # layer
    @abstractmethod
    def plt_section(self, section, x_ticks=None, colors=None, name=None, texts=[]):
        """
        plot the cross section
        :param section: nrow * layer_number or ncol * layer_number
        :param x_ticks: list, e.g. [ 150, 152, 154]
        :param colors: list for colors e.g. ['Sienna', 'SandyBrown', 'SeaGreen', 'Cyan']
        :param name: if not None then save it
        :param texts: [[x, y, text], …]
        :return:
        """
        pass

    @abstractmethod
    def plt_twin(self, _x, _y1, _y2, _xlabel=None, _ylabel1=None, _ylabel2=None, title=None, legend='upper right',
                  name=None, padding_dict={}, texts=[]):
        """
        twin fig
        :param _x: list or array
        :param _y1: list or array
        :param _y2: list or array
        :param _xlabel: str or None
        :param _ylabel1: str or None
        :param _ylabel2: str or None
        :param title: str or None
        :param legend: str
            :parameter upper/lower right/left
        :param name: str if not None then save the fig
        :param padding_dict: left < right  bottom < top        float value adjust the position of the graph
            :parameter left: float default 0.125
            :parameter bottom: float default 0.1
            :parameter right: float default 0.9
            :parameter top: float default 0.9
        :param texts: [[x, y, text], …] or [[x, y, text, color], …]
        :return:
        """
        pass


class Ploter(PlotC):
    def __init__(self):
        super(Ploter, self).__init__()

    @classmethod
    def plt_conterf(cls, matrix, x_ticks=None, y_ticks=None, title=None, colorbar=None, name=None, padding_dict={}):
        nrow, ncol = matrix.shape
        x, y = np.arange(ncol), nrow - np.arange(nrow) - 1
        xx, yy = np.meshgrid(x, y)
        plt.figure()
        plt.subplots_adjust(**padding_dict)
        c = plt.contourf(xx, yy, matrix)
        plt.clabel(c, inline=False, fontsize=10, colors='k')
        if colorbar:
            plt.colorbar(orientation=colorbar)
        if x_ticks is not None:
            plt.xticks([item for item in np.linspace(0, ncol, len(x_ticks))], x_ticks)
        else:
            plt.xticks([])
        if y_ticks is not None:
            plt.yticks([item for item in np.linspace(0, nrow, len(y_ticks))], y_ticks)
        else:
            plt.yticks([])
        if title is not None:
            plt.title(title, fontdict={'fontsize': 16}, backgroundcolor='LightGrey', pad=10)
        if name is None:
            plt.show()
        else:
            plt.savefig(f"{name}.eps", format='eps', dpi=600, bbox_inches='tight')

    @classmethod
    def plt_3d(cls, matrix, x_ticks=None, y_ticks=None, title=None, name=None, padding_dict={}):
        fig = plt.figure()
        ax3 = plt.axes(projection='3d')
        nrow, ncol = matrix.shape
        x, y = np.arange(0, ncol), nrow - np.arange(nrow) - 1
        xx, yy = np.meshgrid(x, y)
        ax3.plot_surface(xx, yy, matrix, cmap='rainbow')
        ax3.auto_scale_xyz([min(x), max(x)], [min(y), max(y)], [matrix.min(), matrix.max()])
        fig.subplots_adjust(**padding_dict)
        # Set the axis labels and titles
        # ax3.set_xlabel(r'$\rho_1$ (g/cm$^3$)')  # Set x lable
        # ax3.set_ylabel(r'$\rho_2$ (g/cm$^3$)')  # Set y lable
        # ax3.set_zlabel('RMSE(m)')               # Set z lable
        if x_ticks is not None:
            plt.xticks([item for item in np.linspace(0, ncol, len(x_ticks))], x_ticks)
        if y_ticks is not None:
            plt.yticks([item for item in np.linspace(0, nrow, len(y_ticks))], y_ticks)
        # ax3.dist = 1.5
        if title is not None:
            plt.title(title, backgroundcolor='LightGrey', pad=10)
        if name is None:
            plt.show()
        else:
            plt.savefig(f"{name}.eps", format='eps', dpi=600)

    @classmethod
    def plt_section(cls, section, x_ticks=None, colors=None, name=None, texts=[]):
        ncol, nchannel = section.shape
        if colors is None:
            colors = ['Sienna', 'SandyBrown', 'Cyan', 'SeaGreen']
        plt.figure()
        x = np.arange(ncol)
        for i in range(nchannel):
            if i == 0:
                layer_up = section[:, i]
                layer_down = np.ones(ncol) * layer_up.min() - 1
                plt.fill_between(x, layer_down, layer_up, color=colors[i], alpha=0.8)
            else:
                layer_down = section[:, i - 1]
                layer_up = section[:, i]
                plt.fill_between(x, layer_down, layer_up, color=colors[i], alpha=0.8)
        if x_ticks is not None:
            plt.xticks([item for item in np.linspace(0, ncol, len(x_ticks))], x_ticks)
        ax = plt.gca()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.ylabel('Elevation(km)')
        for text in texts:
            plt.text(text[0], text[1], text[2])
        if name is None:
            plt.show()
        else:
            plt.savefig(f"{name}.eps", format='eps', dpi=600, bbox_inches='tight')

    @classmethod
    def plt_twin(cls, _x, _y1, _y2, _xlabel=None, _ylabel1=None, _ylabel2=None, title=None, legend='upper right',
                  name=None, padding_dict={}, texts=[]):
        fig = plt.figure()
        plt.subplots_adjust(**padding_dict)
        ax1 = fig.add_subplot()
        color = 'tab:blue'
        if _xlabel:
            ax1.set_xlabel(_xlabel, fontdict={'fontsize': 18})
        if _ylabel1:
            ax1.set_ylabel(_ylabel1, color=color, fontdict={'fontsize': 18})
        line1 = ax1.plot(_x, _y1, color=color, label=_ylabel1, linewidth=2)
        ax1.tick_params(axis='y', labelcolor=color)
        colors = ["tab:blue", "tab:red"] * 10
        for i, text in enumerate(texts):
            plt.text(text[0], text[1], text[2], color=text[3] if len(text) > 3 else colors[i], fontdict={'fontsize': 18})
        ax2 = ax1.twinx()  # 创建共用x轴的第二个y轴

        color = 'tab:red'
        if _ylabel2:
            ax2.set_ylabel(_ylabel2, color=color, fontdict={'fontsize': 18})
        line2 = ax2.plot(_x, _y2, color=color, label=_ylabel2, linestyle='dashed', linewidth=2)
        ax2.tick_params(axis='y', labelcolor=color)

        # fig.tight_layout()
        ax_list = line1 + line2
        labels = [l.get_label() for l in ax_list]
        ax1.legend(ax_list, labels, loc=legend)
        if title is not None:
            plt.title(title, fontdict={'fontsize': 12}, backgroundcolor='LightGrey', pad=10)
        if name is None:
            plt.show()
        else:
            plt.savefig(name)


if __name__ == "__main__":
    x, y = np.linspace(0, 3*np.pi, 101), np.linspace(0, 2*np.pi, 101)
    xx, yy = np.meshgrid(x, y)
    matrix = -4 + 2 * np.sin(xx + yy)
    x_ticks = ["150", "154", "158"]
    y_ticks = ["20", "24", "28"]
    title = "Test"

    ptr = Ploter()

    ptr.plt_conterf(matrix, x_ticks, y_ticks, title, colorbar="vertical")

    ptr.plt_3d(matrix, x_ticks, y_ticks, title)

    section = np.random.uniform(-5, -4, (101, 3))
    section[:, 1] = np.random.uniform(-4, -3, 101)
    section[:, 2] = np.zeros(101)

    ptr.plt_section(section, x_ticks=x_ticks, texts=[[50, -2, "good"]])

    _x = np.linspace(0, 3*np.pi, 100)
    _y1 = np.sin(_x)
    _y2 = 100 * np.cos(_x)

    ptr.plt_twin(_x, _y1, _y2, _xlabel="x", _ylabel1="y1", _ylabel2="y2", title="test",
                 legend="lower left", padding_dict={"right": 0.88}, texts=[[3, 0, "good"], [1, 0.5, "bad"]])
