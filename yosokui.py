import yosoku

import PySimpleGUI as sg

sg.theme('Dark Blue 3')  # please make your windows colorful

layout = [[sg.Text('Enter a Number')],
          [sg.Text('Spacegroup: '), sg.Input('', key='sg')],
          [sg.Text('Unit Cell: '), sg.Input('', key='uc')],
          [sg.Text('High Res: '), sg.Input('', key='hr')],
          [sg.Text('Sequence or # Scatterers: '), sg.Input('', key='ss')],
          [sg.Button('Start'), sg.Button('Exit')]]

window = sg.Window('Yosoku(I) -- Yosoku GUI', layout)

while True:  # Event Loop
    event, values = window.read()
    if event == sg.WIN_CLOSED or event == 'Exit':
        break
    if event == 'Start':
        # change the "output" element to be the value of "input" element
        runpred = yosoku.Phasepred.fromgui(values['sg'], values['uc'], values['hr'], values['ss'])


window.close()

