from PySide2.QtWidgets import *
from main_widget import MainWidget
from main_window import MainWindow

class main:
    def __init__(self):
        pass  # call __init__(self) of the custom base class here

if __name__ == "__main__":

    # Qt Application
    app = QApplication([])

    # main_widget
    app.main_widget = MainWidget()

    # QMainWindow using QWidget as central widget
    app.window = MainWindow(app)

    app.window.show() # shows main window

    app.exec_() # pauses this script until main window is closed