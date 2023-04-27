import bpy
from bpy.types import Operator
from bpy.utils import register_class, unregister_class

class OBJECT_hide_viewport(Operator):
    bl_idname = 'object.hide_viewport'
    bl_label = 'Hide viewport'
    bl_description = 'Globally disable in viewport'
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        for obj in context.selected_objects:
            obj.hide_viewport = True
        return {'FINISHED'}


class OBJECT_hide_viewport_clear(Operator):
    bl_idname = 'object.hide_viewport_clear'
    bl_label = 'Clear viewport hide'
    bl_description = 'Globally cler hiding in viewport'
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        for obj in context.blend_data.objects:
            if obj.hide_viewport:
                obj.hide_viewport = False
        return {'FINISHED'}


def register():
    register_class(OBJECT_hide_viewport)
    register_class(OBJECT_hide_viewport_clear)


def unregister():
    unregister_class(OBJECT_hide_viewport_clear)
    unregister_class(OBJECT_hide_viewport)

if __name__ == '__main__':
    register()