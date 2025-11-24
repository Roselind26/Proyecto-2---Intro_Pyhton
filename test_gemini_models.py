import google.generativeai as genai

# ⚠️ Usa la misma API key que pusiste en el dashboard
API_KEY = "AIzaSyBnvvNgFSEopO1V92teBWrKP-ZZMeZOI3I"

genai.configure(api_key=API_KEY)

try:
    print("Modelos disponibles en tu cuenta:\n")
    for m in genai.list_models():
        # Mostramos nombre y qué métodos soporta
        print(m.name, "->", m.supported_generation_methods)
except Exception as e:
    print("ERROR al llamar a list_models():")
    print(repr(e))
