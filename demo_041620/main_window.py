from PySide2.QtCore import Slot, qApp
from PySide2.QtWidgets import QMainWindow, QAction

class MainWindow(QMainWindow):
    def __init__(self, app):
        QMainWindow.__init__(self)

        # Set up links
        self.app = app
        app.main_window = self
        app.main_widget.app = app
        app.main_widget.main_window = self
        app.main_widget.initialize()

        self.setObjectName('Main Window')
        self.setWindowTitle('Motif Builder Demo 4.16.20')
        
        self.setGeometry(100, 250, 300, 200)
        
        # Exit QAction
        exit_action = QAction("Exit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.exit_app)

        # Window dimensions
        geometry = qApp.desktop().availableGeometry(self)
        self.setFixedSize(geometry.width() * 0.8, geometry.height() * 0.7)
        self.setCentralWidget(app.main_widget)

    @Slot()
    def exit_app(self, checked):
        sys.exit()

    def closeEvent(self, event):

        print('Exiting...')