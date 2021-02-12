import yosoku
import PySimpleGUI as sg
import os


sg.theme("")

inputcolumn = [
    [sg.Text("Enter crystal information and press start")],
    [sg.Text("Spacegroup: "), sg.Input("152", size=(8, 1), key="sg")],
    [sg.Text("Unit Cell: "), sg.Input("158 158 93 90 90 120", key="uc")],
    [sg.Text("High Res: "), sg.Input("2", key="hr")],
    [
        sg.Text(
            """Sequence or 
# Scatterers: """
        ),
        sg.Multiline(
            """DIQMTQTTSSLSASLGDRVTITCSASQGINNYLNWYQQKPDGTVKLLIYYTSSLHSGVPSRFSGSGSGTDYSLTISNLEPEDIATYYCQQYSNLPYTFGGGTKLEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC
        EVQLQQSGPDLVKPGASVKISCKTSGYTFTEYIMHWVKQSHGKSLEWIGGIIPNNGGTSYNQKFKDKATMTVDKSSSTGYMELRSLTSEDSAVYYCTRREVYGRNYYALDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKCCDKTHTCPPCPAPELL
        EPPTACREKQYLINSQCCSLCQPGQKLVSDCTEFTETECLPCGESEFLDTWNRETHCHQHKYCDPNLGLRVQQKGTSETDTICTCEEGWHCTSEACESCVLHRSCSPGFGVKQIATGVSDTICEPCPVGFFSNVSSAFEKCHPWTSCETKDLVVQQAGTNKTDVVCGPQDRLR""",
            size=(50, 10),
            key="ss",
        ),
    ],
    [sg.Button("ASU Predict")],
    [
        sg.Multiline(
            "ASU prediction table will appear here...", size=(60, 8), key="matrup_graph"
        )
    ],
    [
        sg.Text("ASU per mol: "),
        sg.Slider(range=(1, 50), orientation="h", default_value=1, key="am"),
    ],
    [sg.Button("Phase Predict"), sg.Button("Save Plot"), sg.Button("Exit")],
]

plotcolumn = [
    [sg.Image(os.path.join(os.getcwd(), "yosoku.png"), key="predimage")],
]

layout = [[sg.Column(inputcolumn), sg.VSeperator(), sg.Column(plotcolumn),]]

window = sg.Window("Yosoku(I) -- Yosoku GUI", layout)

while True:
    event, values = window.read()
    if event == sg.WIN_CLOSED or event == "Exit":

        break
    if event == "ASU Predict":
        runpred = yosoku.Phasepred(
            values["sg"], values["uc"], values["hr"], values["ss"]
        )
        window["matrup_graph"].update(runpred.asu_pred.table)
        window["am"].update(runpred.asu_pred.n_copies)
    if event == "Phase Predict":
        runpred = yosoku.Phasepred(
            values["sg"], values["uc"], values["hr"], values["ss"]
        )
        runpred.asu_control(values["am"])
        window["predimage"].update(runpred.outfile)
        os.remove(runpred.outfile)

window.close()
