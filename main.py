# main.py
import sys
from PyQt5.QtWidgets import QApplication
from PyQt5.QtCore import Qt  # <--- 导入 Qt 核心命名空间
from main_window import APTKSimulationUI

if __name__ == '__main__':
  
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)

    app = QApplication(sys.argv)

    # 实例化并显示带有顶部主菜单栏的主窗口
    window = APTKSimulationUI()
    window.show()

    # 挂起主事件循环
    sys.exit(app.exec_())